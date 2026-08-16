#!/usr/bin/env python3
"""Hallmark-wide pathway coherence for newly added reverse-consistency results."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from _paths import HALLMARK_GMT, PROJECT_ROOT as ROOT

REVERSE = ROOT / "results" / "extension_reverse_consistency"
OUT = ROOT / "results" / "extension_hallmark_coherence"

SOURCES = {
    "hcl_adult_tissue": ROOT / "results" / "zju_cellatlas_adult_tissue_refined_matching_by_context",
    "hcl_adult_celltype": ROOT / "results" / "zju_cellatlas_adult_celltype_refined_matching_by_context",
}

EDGE_SETS = [
    "trrust",
    "dorothea_ab",
    "trrust_dorothea_intersection",
    "trrust_dorothea_union",
]


def standardize_gene(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper()


def standardize_edges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out["target"] = out["target"].astype(str).str.upper()
    return out


def clean_hallmark_name(name: str) -> str:
    return name.replace("HALLMARK_", "").replace("_", " ").title()


def load_hallmark_gmt(path: Path = HALLMARK_GMT) -> pd.DataFrame:
    rows = []
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, url, *genes = parts
            for gene in genes:
                gene = gene.strip().upper()
                if gene:
                    rows.append(
                        {
                            "collection": "Hallmark",
                            "pathway": clean_hallmark_name(name),
                            "msigdb_name": name,
                            "gene": gene,
                        }
                    )
    df = pd.DataFrame(rows).drop_duplicates()
    if df["msigdb_name"].nunique() != 50:
        raise ValueError(f"Expected 50 Hallmark gene sets, found {df['msigdb_name'].nunique()}")
    return df


def pair_prefix(metric_path: Path) -> str:
    suffix = "_refined_pair_metrics.csv"
    name = metric_path.name
    return name[: -len(suffix)] if name.endswith(suffix) else name


def infer_dataset_condition(pair_prefix_value: str) -> tuple[str, str]:
    if "_celltype_" in pair_prefix_value:
        dataset, condition = pair_prefix_value.split("_celltype_", 1)
        return dataset + "_celltype", condition
    if "_adult_tissue_" in pair_prefix_value:
        dataset, condition = pair_prefix_value.split("_adult_tissue_", 1)
        return dataset + "_adult_tissue", condition
    if "_adult_celltype_" in pair_prefix_value:
        dataset, condition = pair_prefix_value.split("_adult_celltype_", 1)
        return dataset + "_adult_celltype", condition
    parts = pair_prefix_value.split("_")
    return "_".join(parts[:-1]), parts[-1]


def build_universes() -> pd.DataFrame:
    rows = []
    for source_group, source_dir in SOURCES.items():
        for metric_path in sorted(source_dir.rglob("*_refined_pair_metrics.csv")):
            pfx = pair_prefix(metric_path)
            dataset, condition = infer_dataset_condition(pfx)
            for edge_set in EDGE_SETS:
                pos_path = metric_path.with_name(f"{pfx}_{edge_set}_refined_positives.csv")
                neg_path = metric_path.with_name(f"{pfx}_{edge_set}_refined_negatives.csv")
                if not pos_path.exists() or not neg_path.exists():
                    continue
                positives = standardize_edges(pd.read_csv(pos_path, usecols=["tf", "target"]))
                negatives = standardize_edges(pd.read_csv(neg_path, usecols=["repeat", "tf", "target"]))
                positive_targets = set(positives["target"])
                for repeat, group in negatives.groupby("repeat"):
                    universe = positive_targets | set(group["target"])
                    rows.append(
                        {
                            "source_group": source_group,
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


def run_tests(top_pairs: pd.DataFrame, universes: pd.DataFrame, pathways: pd.DataFrame) -> pd.DataFrame:
    pathway_sets = {
        key: set(group["gene"])
        for key, group in pathways.groupby(["collection", "pathway"], sort=False)
    }
    universe_lookup = {
        (row.source_group, row.dataset, row.condition, row.edge_set, int(row.repeat)): set(str(row.universe_targets).split(";"))
        for row in universes.itertuples(index=False)
    }
    rows = []
    group_cols = ["source_group", "dataset", "condition", "edge_set", "repeat", "metric", "top_level"]
    for key, group in top_pairs.groupby(group_cols, sort=False):
        source_group, dataset, condition, edge_set, repeat, metric, top_level = key
        selected = group[group["label"].eq(1)]
        selected_targets = set(standardize_gene(selected["target"]))
        if not selected_targets:
            continue
        universe = universe_lookup[(source_group, dataset, condition, edge_set, int(repeat))]
        universe = {g for g in universe if g and g != "NAN"}
        m = len(universe)
        n = len(selected_targets)
        if m == 0 or n == 0:
            continue
        family = []
        for (collection, pathway), genes in pathway_sets.items():
            geneset = genes & universe
            if not geneset:
                continue
            overlap = selected_targets & geneset
            k = len(overlap)
            expected = n * len(geneset) / m
            fold = (k / n) / (len(geneset) / m) if len(geneset) else np.nan
            pvalue = hypergeom.sf(k - 1, m, len(geneset), n) if k > 0 else 1.0
            family.append(
                {
                    "source_group": source_group,
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "repeat": int(repeat),
                    "metric": metric,
                    "top_level": top_level,
                    "collection": collection,
                    "pathway": pathway,
                    "universe_target_count": m,
                    "selected_target_count": n,
                    "pathway_genes_in_universe": len(geneset),
                    "overlap_genes": k,
                    "expected_overlap": expected,
                    "fold_enrichment": fold,
                    "pvalue": pvalue,
                    "overlap_gene_list": ";".join(sorted(overlap)),
                }
            )
        if family:
            fam = pd.DataFrame(family)
            fam["bh_fdr"] = bh_fdr(fam["pvalue"])
            rows.append(fam)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summarize(tests: pd.DataFrame) -> pd.DataFrame:
    def union_genes(values: pd.Series) -> str:
        genes: set[str] = set()
        for value in values.dropna():
            genes.update(str(value).split(";"))
        return ";".join(sorted(g for g in genes if g))

    return (
        tests.groupby(
            ["source_group", "dataset", "condition", "edge_set", "metric", "top_level", "collection", "pathway"],
            dropna=False,
        )
        .agg(
            repeats_tested=("repeat", "nunique"),
            repeats_with_overlap=("overlap_genes", lambda x: int((x > 0).sum())),
            repeats_fdr_lt_0_05=("bh_fdr", lambda x: int((x < 0.05).sum())),
            mean_selected_targets=("selected_target_count", "mean"),
            median_selected_targets=("selected_target_count", "median"),
            mean_overlap_genes=("overlap_genes", "mean"),
            median_overlap_genes=("overlap_genes", "median"),
            mean_fold_enrichment=("fold_enrichment", "mean"),
            median_fold_enrichment=("fold_enrichment", "median"),
            best_pvalue=("pvalue", "min"),
            best_bh_fdr=("bh_fdr", "min"),
            overlap_gene_union=("overlap_gene_list", union_genes),
        )
        .reset_index()
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pathways = load_hallmark_gmt()
    universes = build_universes()
    top_path = REVERSE / "reverse_extension_top_pairs_by_repeat.csv.gz"
    if not top_path.exists():
        raise FileNotFoundError(top_path)
    top_pairs = pd.read_csv(top_path)
    tests = run_tests(top_pairs, universes, pathways)
    if tests.empty:
        raise RuntimeError("No Hallmark tests were generated.")
    summary = summarize(tests)
    overview = (
        summary.groupby(["source_group", "pathway", "top_level"], dropna=False)
        .agg(
            combinations=("metric", "count"),
            combinations_fdr_in_at_least_5_repeats=("repeats_fdr_lt_0_05", lambda x: int((x >= 5).sum())),
            median_fold_enrichment=("median_fold_enrichment", "median"),
            best_bh_fdr=("best_bh_fdr", "min"),
        )
        .reset_index()
        .sort_values(
            ["source_group", "top_level", "combinations_fdr_in_at_least_5_repeats", "median_fold_enrichment"],
            ascending=[True, True, False, False],
        )
    )
    pathways.to_csv(OUT / "hallmark_full_gene_sets_used.csv", index=False)
    universes.to_csv(OUT / "hallmark_extension_target_universes.csv", index=False)
    tests.to_csv(OUT / "hallmark_extension_tests_by_repeat.csv", index=False)
    summary.to_csv(OUT / "hallmark_extension_summary.csv", index=False)
    overview.to_csv(OUT / "hallmark_extension_overview.csv", index=False)
    report = [
        "# Hallmark full pathway coherence extension report",
        "",
        f"Hallmark GMT: {HALLMARK_GMT}",
        f"Hallmark gene sets: {pathways['pathway'].nunique()}",
        f"Universe rows: {universes.shape[0]}",
        f"Test rows: {tests.shape[0]}",
        "",
        "## Top overview rows",
        "```csv",
        overview.head(100).to_csv(index=False).strip(),
        "```",
    ]
    (OUT / "hallmark_extension_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"hallmark_sets={pathways['pathway'].nunique()}")
    print(f"universe_rows={universes.shape[0]}")
    print(f"test_rows={tests.shape[0]}")
    print(overview.head(60).to_string(index=False))


if __name__ == "__main__":
    main()
