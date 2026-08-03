"""
Reverse validation using already computed refined-pair results.

This is the minimal second objective:
within the TF-gene pairs already evaluated in the strict refined forward
validation, rank pairs by the five existing metrics and ask how many top-ranked
pairs are curated positives.

No new pair universe is generated and no expression metrics are recomputed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[1]
FORWARD = ROOT / "results" / "forward_validation_refined"
LOCAL_REFINED = ROOT / "results" / "formal_abs_association_refined_matching"
SERVER_REFINED = ROOT / "results" / "server_added_refined_matching"
OUT = ROOT / "results" / "reverse_validation_from_refined_pairs"

METRICS = [
    "pearson",
    "spearman",
    "mutual_information",
    "coexpression_probability",
    "codetection_odds_ratio",
]
TOP_LEVELS = {
    "top_100": 100,
    "top_500": 500,
    "top_5pct": "5pct",
}


SOURCE_DIRS = {
    "local_refined": LOCAL_REFINED,
    "server_nygc": SERVER_REFINED / "nygc_multimodal_pbmc_refined_matching",
    "server_gse126030": SERVER_REFINED / "gse126030_refined_matching",
}


def standardize_edges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out["target"] = out["target"].astype(str).str.upper()
    return out


def load_known_sets() -> dict[str, set[tuple[str, str]]]:
    trrust = standardize_edges(pd.read_csv(ROOT / "data" / "ground_truth" / "trrust_edges_standardized.csv"))
    dorothea = standardize_edges(pd.read_csv(ROOT / "data" / "ground_truth" / "dorothea_ab_edges_standardized.csv"))
    trrust_set = set(map(tuple, trrust[["tf", "target"]].drop_duplicates().itertuples(index=False, name=None)))
    dorothea_set = set(map(tuple, dorothea[["tf", "target"]].drop_duplicates().itertuples(index=False, name=None)))
    return {
        "trrust": trrust_set,
        "dorothea_ab": dorothea_set,
        "trrust_dorothea_intersection": trrust_set & dorothea_set,
        "trrust_dorothea_union": trrust_set | dorothea_set,
    }


def load_pathways() -> pd.DataFrame:
    paths = [
        ROOT / "data" / "ground_truth" / "pathway_gene_sets" / "kang_ifn_relevant_pathways.csv",
        ROOT / "data" / "ground_truth" / "pathway_gene_sets" / "kang_reactome_ifn_relevant_pathways.csv",
    ]
    frames = []
    for path in paths:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df["gene"] = df["gene"].astype(str).str.upper()
        df["collection"] = "Reactome" if "reactome" in path.name.lower() else "Hallmark"
        frames.append(df[["collection", "pathway", "gene"]].drop_duplicates())
    return pd.concat(frames, ignore_index=True).drop_duplicates() if frames else pd.DataFrame()


def source_dir_for(result_source: str) -> Path:
    if result_source not in SOURCE_DIRS:
        raise KeyError(f"Unknown result_source: {result_source}")
    return SOURCE_DIRS[result_source]


def prefix(dataset: str, condition: str, edge_set: str | None = None) -> str:
    parts = [dataset, condition]
    if edge_set is not None:
        parts.append(edge_set)
    return "_".join(parts).replace("/", "_")


def annotate_known(df: pd.DataFrame, known_sets: dict[str, set[tuple[str, str]]]) -> pd.DataFrame:
    out = df.copy()
    pairs = list(map(tuple, out[["tf", "target"]].itertuples(index=False, name=None)))
    for name, known in known_sets.items():
        out[f"in_{name}"] = [p in known for p in pairs]
    return out


def load_unit_pairs(
    source_dir: Path,
    dataset: str,
    condition: str,
    edge_set: str,
    known_sets: dict[str, set[tuple[str, str]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    unit_prefix = prefix(dataset, condition, edge_set)
    pair_prefix = prefix(dataset, condition)
    pos_path = source_dir / f"{unit_prefix}_refined_positives.csv"
    neg_path = source_dir / f"{unit_prefix}_refined_negatives.csv"
    metric_path = source_dir / f"{pair_prefix}_refined_pair_metrics.csv"
    for path in [pos_path, neg_path, metric_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    positives = standardize_edges(pd.read_csv(pos_path))
    negatives = standardize_edges(pd.read_csv(neg_path))
    metrics = standardize_edges(pd.read_csv(metric_path))
    metrics = metrics.drop_duplicates(["tf", "target"]).reset_index(drop=True)
    for metric in ["pearson", "spearman"]:
        metrics[metric] = metrics[metric].abs()

    positives = positives[["tf", "target"]].drop_duplicates().copy()
    positives["label"] = 1
    positives["repeat"] = -1
    positives["edge_role"] = "known_positive"
    negatives = negatives[["repeat", "tf", "target"]].drop_duplicates().copy()
    negatives["label"] = 0
    negatives["edge_role"] = "matched_background"

    all_pairs = pd.concat([positives, negatives], ignore_index=True)
    all_pairs = all_pairs.merge(metrics, on=["tf", "target"], how="left")
    all_pairs = annotate_known(all_pairs, known_sets)
    return all_pairs, metrics


def top_k_for_level(total: int, level: str, value: int | str) -> int:
    if value == "5pct":
        return max(1, int(np.ceil(total * 0.05)))
    return min(int(value), total)


def analyze_unit(
    result_source: str,
    dataset: str,
    condition: str,
    edge_set: str,
    known_sets: dict[str, set[tuple[str, str]]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source_dir = source_dir_for(result_source)
    pairs, _metrics = load_unit_pairs(source_dir, dataset, condition, edge_set, known_sets)
    positives = pairs[pairs["label"].eq(1)].copy()
    negatives = pairs[pairs["label"].eq(0)].copy()
    repeats = sorted(negatives["repeat"].dropna().astype(int).unique())

    search_rows = []
    top_rows = []
    overlap_rows = []

    baseline_known = {}
    positive_pairs = set(map(tuple, positives[["tf", "target"]].itertuples(index=False, name=None)))
    negative_pairs = set(map(tuple, negatives[["tf", "target"]].drop_duplicates().itertuples(index=False, name=None)))
    pair_universe = positive_pairs | negative_pairs
    for known_name, known in known_sets.items():
        baseline_known[known_name] = len(pair_universe & known)

    search_rows.append(
        {
            "result_source": result_source,
            "dataset": dataset,
            "condition": condition,
            "edge_set": edge_set,
            "positive_pairs": len(positive_pairs),
            "unique_negative_pairs_across_repeats": len(negative_pairs),
            "unique_pair_universe": len(pair_universe),
            "per_repeat_expected_positive_fraction": 0.5,
            **{f"{name}_known_pairs_in_universe": count for name, count in baseline_known.items()},
            **{f"{name}_known_fraction_in_universe": count / len(pair_universe) for name, count in baseline_known.items()},
        }
    )

    for repeat in repeats:
        repeat_neg = negatives[negatives["repeat"].eq(repeat)]
        test = pd.concat(
            [
                positives.assign(repeat=repeat),
                repeat_neg,
            ],
            ignore_index=True,
        )
        for metric in METRICS:
            valid = test.replace([np.inf, -np.inf], np.nan).dropna(subset=[metric]).copy()
            if valid.empty:
                continue
            valid = valid.sort_values(metric, ascending=False).reset_index(drop=True)
            valid["rank"] = np.arange(1, len(valid) + 1)
            total = len(valid)
            for level, raw_k in TOP_LEVELS.items():
                k = top_k_for_level(total, level, raw_k)
                sub = valid.head(k).copy()
                sub["result_source"] = result_source
                sub["dataset"] = dataset
                sub["condition"] = condition
                sub["edge_set"] = edge_set
                sub["metric"] = metric
                sub["top_level"] = level
                sub["top_n"] = k
                cols = [
                    "result_source",
                    "dataset",
                    "condition",
                    "edge_set",
                    "repeat",
                    "metric",
                    "top_level",
                    "top_n",
                    "rank",
                    "tf",
                    "target",
                    "label",
                    "edge_role",
                    "pearson",
                    "spearman",
                    "mutual_information",
                    "coexpression_probability",
                    "codetection_odds_ratio",
                    "in_trrust",
                    "in_dorothea_ab",
                    "in_trrust_dorothea_intersection",
                    "in_trrust_dorothea_union",
                ]
                top_rows.append(sub[cols])
                overlap_rows.append(
                    {
                        "result_source": result_source,
                        "dataset": dataset,
                        "condition": condition,
                        "edge_set": edge_set,
                        "repeat": repeat,
                        "metric": metric,
                        "top_level": level,
                        "top_n": k,
                        "known_positive_in_top": int(sub["label"].sum()),
                        "known_positive_fraction_in_top": float(sub["label"].mean()),
                        "balanced_random_baseline": 0.5,
                        "delta_vs_balanced_baseline": float(sub["label"].mean() - 0.5),
                        "trrust_union_known_in_top": int(sub["in_trrust_dorothea_union"].sum()),
                        "trrust_union_known_fraction_in_top": float(sub["in_trrust_dorothea_union"].mean()),
                    }
                )

    return pd.DataFrame(search_rows), pd.DataFrame(overlap_rows), pd.concat(top_rows, ignore_index=True)


def summarize_overlap(overlap: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        overlap.groupby(["result_source", "dataset", "condition", "edge_set", "metric", "top_level"], dropna=False)
        .agg(
            n_repeats=("repeat", "nunique"),
            top_n_mean=("top_n", "mean"),
            known_positive_fraction_mean=("known_positive_fraction_in_top", "mean"),
            known_positive_fraction_sd=("known_positive_fraction_in_top", "std"),
            delta_vs_balanced_baseline_mean=("delta_vs_balanced_baseline", "mean"),
            repeats_above_0_5=("known_positive_fraction_in_top", lambda x: int((x > 0.5).sum())),
            trrust_union_known_fraction_mean=("trrust_union_known_fraction_in_top", "mean"),
        )
        .reset_index()
    )
    return grouped.sort_values(["dataset", "condition", "edge_set", "metric", "top_level"]).reset_index(drop=True)


def consensus_pairs(top_pairs: pd.DataFrame) -> pd.DataFrame:
    top500 = top_pairs[top_pairs["top_level"].eq("top_500")].copy()
    if top500.empty:
        return pd.DataFrame()
    metrics = METRICS
    score_agg = {m: (m, "mean") for m in metrics}
    out = (
        top500.groupby(["dataset", "condition", "edge_set", "metric", "tf", "target", "label"], dropna=False)
        .agg(
            repeat_top_count=("repeat", "nunique"),
            best_rank=("rank", "min"),
            mean_rank=("rank", "mean"),
            **score_agg,
            in_trrust=("in_trrust", "max"),
            in_dorothea_ab=("in_dorothea_ab", "max"),
            in_trrust_dorothea_intersection=("in_trrust_dorothea_intersection", "max"),
            in_trrust_dorothea_union=("in_trrust_dorothea_union", "max"),
        )
        .reset_index()
    )
    return out.sort_values(
        ["dataset", "condition", "edge_set", "metric", "repeat_top_count", "best_rank"],
        ascending=[True, True, True, True, False, True],
    ).reset_index(drop=True)


def pathway_summary(top_pairs: pd.DataFrame, pathways: pd.DataFrame) -> pd.DataFrame:
    if pathways.empty or top_pairs.empty:
        return pd.DataFrame()
    rows = []
    for (dataset, condition, edge_set, metric, top_level), group in top_pairs.groupby(
        ["dataset", "condition", "edge_set", "metric", "top_level"]
    ):
        if top_level not in {"top_100", "top_500"}:
            continue
        targets = set(group["target"].astype(str).str.upper())
        universe = set(top_pairs[(top_pairs["dataset"].eq(dataset)) & (top_pairs["condition"].eq(condition))]["target"].astype(str).str.upper())
        m = len(universe)
        n = len(targets)
        if m == 0 or n == 0:
            continue
        for (collection, pathway), pgroup in pathways.groupby(["collection", "pathway"]):
            geneset = set(pgroup["gene"]) & universe
            if not geneset:
                continue
            overlap = targets & geneset
            if not overlap:
                continue
            k = len(overlap)
            pvalue = hypergeom.sf(k - 1, m, len(geneset), n)
            rows.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "metric": metric,
                    "top_level": top_level,
                    "collection": collection,
                    "pathway": pathway,
                    "top_target_genes": n,
                    "pathway_genes_in_universe": len(geneset),
                    "overlap_genes": k,
                    "fold_enrichment": (k / n) / (len(geneset) / m),
                    "pvalue": pvalue,
                    "overlap_gene_list": ";".join(sorted(overlap)),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values("pvalue").reset_index(drop=True)
    out["bh_fdr"] = np.minimum.accumulate((out["pvalue"] * len(out) / (np.arange(len(out)) + 1))[::-1])[::-1]
    return out.sort_values(["dataset", "condition", "edge_set", "metric", "top_level", "bh_fdr"]).reset_index(drop=True)


def write_report(search: pd.DataFrame, overlap_summary: pd.DataFrame, top_pairs: pd.DataFrame, pathways: pd.DataFrame) -> None:
    focus = overlap_summary[
        overlap_summary["metric"].isin(["pearson", "spearman"])
        & overlap_summary["top_level"].isin(["top_100", "top_500"])
    ].copy()
    lines = [
        "# Reverse Validation From Refined Pairs",
        "",
        "Scope: rank only the TF-gene pairs already evaluated in the strict refined forward validation.",
        "No new candidate universe is generated and no expression metric is recomputed.",
        "",
        f"- Analysis units: {search[['dataset', 'condition', 'edge_set']].drop_duplicates().shape[0]}",
        f"- Top-pair rows: {len(top_pairs)}",
        f"- Overlap summary rows: {len(overlap_summary)}",
        f"- Pathway summary rows: {len(pathways)}",
        "",
        "Balanced random baseline: each repeat contains 1:1 curated positives and matched background pairs, so the expected curated-positive fraction among top-ranked pairs is 0.5 if the metric carries no ranking signal.",
        "",
        "Pearson/Spearman top curated-positive fractions:",
        focus[
            [
                "dataset",
                "condition",
                "edge_set",
                "metric",
                "top_level",
                "known_positive_fraction_mean",
                "delta_vs_balanced_baseline_mean",
                "repeats_above_0_5",
            ]
        ].to_string(index=False),
    ]
    (OUT / "reverse_refined_validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    full_summary = pd.read_csv(FORWARD / "forward_refined_full_summary.csv")
    units = full_summary[["result_source", "dataset", "condition", "edge_set"]].drop_duplicates().reset_index(drop=True)
    known_sets = load_known_sets()
    pathways = load_pathways()

    search_rows = []
    overlap_rows = []
    top_rows = []
    for _, row in units.iterrows():
        print(f"[unit] {row['dataset']} {row['condition']} {row['edge_set']}", flush=True)
        search, overlap, top = analyze_unit(
            row["result_source"],
            row["dataset"],
            row["condition"],
            row["edge_set"],
            known_sets,
        )
        search_rows.append(search)
        overlap_rows.append(overlap)
        top_rows.append(top)

    search_all = pd.concat(search_rows, ignore_index=True)
    overlap_by_repeat = pd.concat(overlap_rows, ignore_index=True)
    top_all = pd.concat(top_rows, ignore_index=True)
    overlap_summary = summarize_overlap(overlap_by_repeat)
    consensus = consensus_pairs(top_all)
    pathway = pathway_summary(top_all, pathways)

    search_all.to_csv(OUT / "reverse_refined_search_space_summary.csv", index=False)
    overlap_by_repeat.to_csv(OUT / "reverse_refined_known_overlap_by_repeat.csv", index=False)
    overlap_summary.to_csv(OUT / "reverse_refined_known_overlap_summary.csv", index=False)
    top_all.to_csv(OUT / "reverse_refined_top_pairs_by_repeat.csv", index=False)
    consensus.to_csv(OUT / "reverse_refined_consensus_top_pairs.csv", index=False)
    pathway.to_csv(OUT / "reverse_refined_pathway_summary.csv", index=False)
    write_report(search_all, overlap_summary, top_all, pathway)
    print(f"[write] {OUT}", flush=True)


if __name__ == "__main__":
    main()
