"""
Pathway validation for the simplified reverse analysis.

This uses only high-scoring TF-gene pairs already selected from the refined
positive/background pair space. For each matching repeat, top target genes are
tested against the target-gene universe available in that same repeat.
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
REVERSE = ROOT / "results" / "reverse_validation_from_refined_pairs"

SOURCE_DIRS = {
    "local_refined": LOCAL_REFINED,
    "server_nygc": SERVER_REFINED / "nygc_multimodal_pbmc_refined_matching",
    "server_gse126030": SERVER_REFINED / "gse126030_refined_matching",
}

TOP_LEVELS = ["top_100", "top_500", "top_5pct"]


def prefix(dataset: str, condition: str, edge_set: str | None = None) -> str:
    parts = [dataset, condition]
    if edge_set is not None:
        parts.append(edge_set)
    return "_".join(parts).replace("/", "_")


def standardize_gene_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper()


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
        if not {"pathway", "gene"}.issubset(df.columns):
            raise ValueError(f"Pathway file missing required columns: {path}")
        df = df[["pathway", "gene"]].copy()
        df["gene"] = standardize_gene_series(df["gene"])
        df["collection"] = "Reactome" if "reactome" in path.name.lower() else "Hallmark"
        frames.append(df[["collection", "pathway", "gene"]])
    if not frames:
        return pd.DataFrame(columns=["collection", "pathway", "gene"])
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def load_repeat_universes() -> pd.DataFrame:
    full = pd.read_csv(FORWARD / "forward_refined_full_summary.csv")
    units = full[["result_source", "dataset", "condition", "edge_set"]].drop_duplicates()
    rows = []
    for _, unit in units.iterrows():
        result_source = unit["result_source"]
        source_dir = SOURCE_DIRS[result_source]
        dataset = unit["dataset"]
        condition = unit["condition"]
        edge_set = unit["edge_set"]
        unit_prefix = prefix(dataset, condition, edge_set)
        pos_path = source_dir / f"{unit_prefix}_refined_positives.csv"
        neg_path = source_dir / f"{unit_prefix}_refined_negatives.csv"
        if not pos_path.exists() or not neg_path.exists():
            raise FileNotFoundError(f"Missing refined pairs for {unit_prefix}")
        positives = pd.read_csv(pos_path, usecols=["target"])
        negatives = pd.read_csv(neg_path, usecols=["repeat", "target"])
        positive_targets = set(standardize_gene_series(positives["target"]))
        negatives["target"] = standardize_gene_series(negatives["target"])
        for repeat, group in negatives.groupby("repeat"):
            universe = positive_targets | set(group["target"])
            rows.append(
                {
                    "result_source": result_source,
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "repeat": int(repeat),
                    "universe_targets": ";".join(sorted(universe)),
                    "universe_target_count": len(universe),
                    "positive_target_count": len(positive_targets),
                    "negative_target_count": group["target"].nunique(),
                }
            )
    return pd.DataFrame(rows)


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pvalues.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    q_sorted = ranked * m / (np.arange(m) + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return pd.Series(q, index=pvalues.index)


def pathway_tests(top_pairs: pd.DataFrame, universes: pd.DataFrame, pathways: pd.DataFrame) -> pd.DataFrame:
    pathway_groups = {
        key: set(group["gene"])
        for key, group in pathways.groupby(["collection", "pathway"], sort=False)
    }
    universe_lookup = {
        (
            row.dataset,
            row.condition,
            row.edge_set,
            int(row.repeat),
        ): set(str(row.universe_targets).split(";"))
        for row in universes.itertuples(index=False)
    }
    rows = []
    group_cols = ["dataset", "condition", "edge_set", "repeat", "metric", "top_level"]
    for key, group in top_pairs.groupby(group_cols, sort=False):
        dataset, condition, edge_set, repeat, metric, top_level = key
        universe = universe_lookup[(dataset, condition, edge_set, int(repeat))]
        top_targets = set(standardize_gene_series(group["target"]))
        m = len(universe)
        n = len(top_targets)
        if m == 0 or n == 0:
            continue
        for (collection, pathway), genes in pathway_groups.items():
            geneset = genes & universe
            if not geneset:
                continue
            overlap = top_targets & geneset
            k = len(overlap)
            expected = n * len(geneset) / m
            fold = (k / n) / (len(geneset) / m) if len(geneset) else np.nan
            pvalue = hypergeom.sf(k - 1, m, len(geneset), n) if k > 0 else 1.0
            rows.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "repeat": int(repeat),
                    "metric": metric,
                    "top_level": top_level,
                    "collection": collection,
                    "pathway": pathway,
                    "universe_target_count": m,
                    "top_target_count": n,
                    "pathway_genes_in_universe": len(geneset),
                    "overlap_genes": k,
                    "expected_overlap": expected,
                    "fold_enrichment": fold,
                    "pvalue": pvalue,
                    "overlap_gene_list": ";".join(sorted(overlap)),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["bh_fdr"] = bh_fdr(out["pvalue"])
    return out.sort_values(
        ["dataset", "condition", "edge_set", "metric", "top_level", "bh_fdr", "pvalue"]
    ).reset_index(drop=True)


def summarize_tests(tests: pd.DataFrame) -> pd.DataFrame:
    if tests.empty:
        return tests

    def union_genes(values: pd.Series) -> str:
        genes: set[str] = set()
        for value in values.dropna():
            if value:
                genes.update(str(value).split(";"))
        return ";".join(sorted(g for g in genes if g))

    grouped = (
        tests.groupby(
            ["dataset", "condition", "edge_set", "metric", "top_level", "collection", "pathway"],
            dropna=False,
        )
        .agg(
            repeats_tested=("repeat", "nunique"),
            repeats_with_overlap=("overlap_genes", lambda x: int((x > 0).sum())),
            repeats_fdr_lt_0_05=("bh_fdr", lambda x: int((x < 0.05).sum())),
            mean_top_targets=("top_target_count", "mean"),
            mean_pathway_genes_in_universe=("pathway_genes_in_universe", "mean"),
            mean_overlap_genes=("overlap_genes", "mean"),
            median_overlap_genes=("overlap_genes", "median"),
            mean_expected_overlap=("expected_overlap", "mean"),
            mean_fold_enrichment=("fold_enrichment", "mean"),
            median_fold_enrichment=("fold_enrichment", "median"),
            best_pvalue=("pvalue", "min"),
            best_bh_fdr=("bh_fdr", "min"),
            overlap_gene_union=("overlap_gene_list", union_genes),
        )
        .reset_index()
    )
    return grouped.sort_values(
        [
            "dataset",
            "condition",
            "edge_set",
            "metric",
            "top_level",
            "repeats_fdr_lt_0_05",
            "mean_fold_enrichment",
        ],
        ascending=[True, True, True, True, True, False, False],
    ).reset_index(drop=True)


def write_report(summary: pd.DataFrame, tests: pd.DataFrame, universes: pd.DataFrame) -> None:
    focus = summary[
        summary["metric"].isin(["pearson", "spearman"])
        & summary["top_level"].isin(["top_100", "top_500"])
    ].copy()
    top = focus.sort_values(
        ["repeats_fdr_lt_0_05", "mean_fold_enrichment", "mean_overlap_genes"],
        ascending=[False, False, False],
    ).head(30)
    lines = [
        "# Reverse Pathway Validation From Refined Pairs",
        "",
        "Scope: pathway validation is performed only on high-scoring pairs from the existing refined positive/background pair space.",
        "The universe for each test is the target-gene set available in the same dataset, condition, edge set, and matching repeat.",
        "",
        f"- Repeat-level pathway tests: {len(tests)}",
        f"- Summary rows: {len(summary)}",
        f"- Analysis units with repeat universes: {universes[['dataset', 'condition', 'edge_set']].drop_duplicates().shape[0]}",
        "",
        "Interpretation: pathway enrichment supports functional consistency of high-scoring targets; it does not prove causal regulation.",
        "",
        "Top Pearson/Spearman pathway signals:",
        top[
            [
                "dataset",
                "condition",
                "edge_set",
                "metric",
                "top_level",
                "collection",
                "pathway",
                "repeats_with_overlap",
                "repeats_fdr_lt_0_05",
                "mean_overlap_genes",
                "mean_fold_enrichment",
                "best_bh_fdr",
            ]
        ].to_string(index=False),
    ]
    (REVERSE / "reverse_refined_pathway_validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    REVERSE.mkdir(parents=True, exist_ok=True)
    top_pair_path = REVERSE / "reverse_refined_top_pairs_by_repeat.csv"
    if not top_pair_path.exists():
        required = [
            REVERSE / "reverse_refined_pathway_validation_summary.csv",
            REVERSE / "reverse_refined_pathway_validation_by_repeat.csv",
            REVERSE / "reverse_refined_pathway_repeat_universes.csv",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing pathway validation inputs. Either restore "
                "reverse_refined_top_pairs_by_repeat.csv or provide existing "
                f"pathway summary outputs. Missing files: {missing}"
            )
        summary = pd.read_csv(required[0])
        tests = pd.read_csv(required[1])
        universes = pd.read_csv(required[2])
        print("[check] reverse_refined_top_pairs_by_repeat.csv is not included in this GitHub package.")
        print("[check] Existing pathway validation outputs are present and readable.")
        print(f"[summary rows] {len(summary)}")
        print(f"[repeat-level test rows] {len(tests)}")
        print(f"[repeat universe rows] {len(universes)}")
        return

    pathways = load_pathways()
    universes = load_repeat_universes()
    usecols = ["dataset", "condition", "edge_set", "repeat", "metric", "top_level", "target"]
    top_pairs = pd.read_csv(top_pair_path, usecols=usecols)
    top_pairs = top_pairs[top_pairs["top_level"].isin(TOP_LEVELS)].copy()
    top_pairs["target"] = standardize_gene_series(top_pairs["target"])
    tests = pathway_tests(top_pairs, universes, pathways)
    summary = summarize_tests(tests)
    universes.to_csv(REVERSE / "reverse_refined_pathway_repeat_universes.csv", index=False)
    tests.to_csv(REVERSE / "reverse_refined_pathway_validation_by_repeat.csv", index=False)
    summary.to_csv(REVERSE / "reverse_refined_pathway_validation_summary.csv", index=False)
    write_report(summary, tests, universes)
    print(f"[write] {REVERSE}")


if __name__ == "__main__":
    main()
