#!/usr/bin/env python3
"""Reverse consistency for newly added refined-matching contexts.

This script reuses already computed refined positive/background pairs and
pair-level metric scores. It does not recompute expression metrics or generate
any new pair universe.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from _paths import PROJECT_ROOT as ROOT, project_relative

OUT = ROOT / "results" / "extension_reverse_consistency"

SOURCES = {
    "hcl_adult_tissue": ROOT / "results" / "zju_cellatlas_adult_tissue_refined_matching_by_context",
    "hcl_adult_celltype": ROOT / "results" / "zju_cellatlas_adult_celltype_refined_matching_by_context",
}

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

EDGE_SETS = [
    "trrust",
    "dorothea_ab",
    "trrust_dorothea_intersection",
    "trrust_dorothea_union",
]


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)


def standardize_edges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out["target"] = out["target"].astype(str).str.upper()
    return out


def discover_metric_files(source_dir: Path) -> list[Path]:
    return sorted(source_dir.rglob("*_refined_pair_metrics.csv"))


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


def load_unit(metric_path: Path, edge_set: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    pfx = pair_prefix(metric_path)
    pos_path = metric_path.with_name(f"{pfx}_{edge_set}_refined_positives.csv")
    neg_path = metric_path.with_name(f"{pfx}_{edge_set}_refined_negatives.csv")
    if not pos_path.exists() or not neg_path.exists():
        return pd.DataFrame(), pd.DataFrame()

    metrics = standardize_edges(read_csv(metric_path))
    metrics = metrics.drop_duplicates(["tf", "target"])
    for metric in ["pearson", "spearman"]:
        if metric in metrics.columns:
            metrics[metric] = metrics[metric].abs()

    positives = standardize_edges(read_csv(pos_path))
    positives = positives[["tf", "target"]].drop_duplicates()
    positives["label"] = 1
    positives["edge_role"] = "curated_positive"

    negatives = standardize_edges(read_csv(neg_path))
    negatives = negatives[["repeat", "tf", "target"]].drop_duplicates()
    negatives["label"] = 0
    negatives["edge_role"] = "matched_background"

    positives = positives.merge(metrics, on=["tf", "target"], how="left")
    negatives = negatives.merge(metrics, on=["tf", "target"], how="left")
    return positives, negatives


def top_count(total_pairs: int, level: str, raw_value: int | str) -> tuple[int, bool]:
    if raw_value == "5pct":
        return max(1, int(np.ceil(total_pairs * 0.05))), True
    requested = int(raw_value)
    return min(requested, total_pairs), total_pairs >= requested


def sort_ranked_pairs(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    return df.sort_values(metric, ascending=False, kind="quicksort").reset_index(drop=True)


def ci95(values: pd.Series) -> float:
    vals = values.dropna().to_numpy(dtype=float)
    if len(vals) <= 1:
        return np.nan
    return 1.96 * vals.std(ddof=1) / np.sqrt(len(vals))


def analyze() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    by_repeat_rows: list[dict] = []
    top_pair_chunks: list[pd.DataFrame] = []
    unit_inventory_rows: list[dict] = []

    for source_group, source_dir in SOURCES.items():
        metric_files = discover_metric_files(source_dir)
        for metric_path in metric_files:
            pfx = pair_prefix(metric_path)
            dataset, condition = infer_dataset_condition(pfx)
            for edge_set in EDGE_SETS:
                positives, negatives = load_unit(metric_path, edge_set)
                if positives.empty or negatives.empty:
                    continue
                repeats = sorted(negatives["repeat"].dropna().astype(int).unique())
                unit_inventory_rows.append(
                    {
                        "source_group": source_group,
                        "dataset": dataset,
                        "condition": condition,
                        "edge_set": edge_set,
                        "positive_pairs": positives[["tf", "target"]].drop_duplicates().shape[0],
                        "unique_negative_pairs_across_repeats": negatives[["tf", "target"]].drop_duplicates().shape[0],
                        "repeats": len(repeats),
                        "metric_file": project_relative(metric_path),
                    }
                )
                for repeat in repeats:
                    repeat_neg = negatives[negatives["repeat"].astype(int).eq(int(repeat))].copy()
                    test = pd.concat([positives.assign(repeat=int(repeat)), repeat_neg], ignore_index=True)
                    pos_n = int(test["label"].sum())
                    neg_n = int((test["label"] == 0).sum())
                    for metric in METRICS:
                        valid = test.replace([np.inf, -np.inf], np.nan).dropna(subset=[metric]).copy()
                        if valid.empty:
                            continue
                        valid = sort_ranked_pairs(valid, metric)
                        valid["rank"] = np.arange(1, len(valid) + 1)
                        total_pairs = len(valid)
                        for top_level, raw_top in TOP_LEVELS.items():
                            actual_top_n, complete_threshold = top_count(total_pairs, top_level, raw_top)
                            sub = valid.head(actual_top_n).copy()
                            pos_top = int(sub["label"].sum())
                            frac = float(sub["label"].mean())
                            by_repeat_rows.append(
                                {
                                    "source_group": source_group,
                                    "dataset": dataset,
                                    "condition": condition,
                                    "edge_set": edge_set,
                                    "repeat": int(repeat),
                                    "metric": metric,
                                    "top_level": top_level,
                                    "requested_top_n": raw_top,
                                    "actual_top_n": int(actual_top_n),
                                    "complete_threshold": bool(complete_threshold),
                                    "total_pairs_ranked": int(total_pairs),
                                    "positive_pairs_in_evaluation": pos_n,
                                    "negative_pairs_in_evaluation": neg_n,
                                    "curated_positive_in_top": pos_top,
                                    "curated_positive_fraction": frac,
                                    "balanced_random_baseline": 0.5,
                                    "delta_vs_0_5": frac - 0.5,
                                }
                            )
                            export_cols = [
                                "repeat",
                                "rank",
                                "tf",
                                "target",
                                "label",
                                "edge_role",
                                metric,
                            ]
                            tmp = sub[export_cols].copy()
                            tmp["source_group"] = source_group
                            tmp["dataset"] = dataset
                            tmp["condition"] = condition
                            tmp["edge_set"] = edge_set
                            tmp["metric"] = metric
                            tmp["top_level"] = top_level
                            tmp["actual_top_n"] = int(actual_top_n)
                            top_pair_chunks.append(tmp)

    inventory = pd.DataFrame(unit_inventory_rows)
    by_repeat = pd.DataFrame(by_repeat_rows)
    if by_repeat.empty:
        raise RuntimeError("No reverse-consistency rows were generated.")

    summary = (
        by_repeat.groupby(["source_group", "dataset", "condition", "edge_set", "metric", "top_level"], dropna=False)
        .agg(
            repeats=("repeat", "nunique"),
            mean_curated_positive_fraction=("curated_positive_fraction", "mean"),
            sd_curated_positive_fraction=("curated_positive_fraction", "std"),
            ci95_curated_positive_fraction=("curated_positive_fraction", ci95),
            min_curated_positive_fraction=("curated_positive_fraction", "min"),
            max_curated_positive_fraction=("curated_positive_fraction", "max"),
            repeats_above_0_5=("curated_positive_fraction", lambda x: int((x > 0.5).sum())),
            mean_delta_vs_0_5=("delta_vs_0_5", "mean"),
            mean_actual_top_n=("actual_top_n", "mean"),
            complete_threshold_repeats=("complete_threshold", lambda x: int(pd.Series(x).sum())),
            mean_total_pairs_ranked=("total_pairs_ranked", "mean"),
            mean_positive_pairs_in_evaluation=("positive_pairs_in_evaluation", "mean"),
            mean_negative_pairs_in_evaluation=("negative_pairs_in_evaluation", "mean"),
        )
        .reset_index()
    )
    summary["all_repeats_above_0_5"] = summary["repeats_above_0_5"].eq(summary["repeats"])

    overview = (
        summary.groupby(["source_group", "metric", "top_level"], dropna=False)
        .agg(
            analysis_units=("condition", "count"),
            median_unit_mean_fraction=("mean_curated_positive_fraction", "median"),
            mean_unit_mean_fraction=("mean_curated_positive_fraction", "mean"),
            min_unit_mean_fraction=("mean_curated_positive_fraction", "min"),
            max_unit_mean_fraction=("mean_curated_positive_fraction", "max"),
            units_mean_above_0_5=("mean_curated_positive_fraction", lambda x: int((x > 0.5).sum())),
            units_all_repeats_above_0_5=("all_repeats_above_0_5", lambda x: int(pd.Series(x).sum())),
            units_with_complete_threshold=("complete_threshold_repeats", lambda x: int((x == 10).sum())),
            median_actual_top_n=("mean_actual_top_n", "median"),
        )
        .reset_index()
    )

    by_repeat.to_csv(OUT / "reverse_extension_known_fraction_by_repeat.csv", index=False)
    summary.to_csv(OUT / "reverse_extension_known_fraction_summary.csv", index=False)
    overview.to_csv(OUT / "reverse_extension_metric_overview.csv", index=False)
    inventory.to_csv(OUT / "reverse_extension_unit_inventory.csv", index=False)

    top_pairs = pd.concat(top_pair_chunks, ignore_index=True) if top_pair_chunks else pd.DataFrame()
    if not top_pairs.empty:
        top_pairs.to_csv(OUT / "reverse_extension_top_pairs_by_repeat.csv.gz", index=False, compression="gzip")
        consensus = (
            top_pairs.groupby(["source_group", "dataset", "condition", "edge_set", "metric", "top_level", "tf", "target", "label"], dropna=False)
            .agg(
                times_selected=("repeat", "nunique"),
                best_rank=("rank", "min"),
                median_rank=("rank", "median"),
            )
            .reset_index()
            .sort_values(["source_group", "dataset", "condition", "edge_set", "metric", "top_level", "best_rank", "median_rank"])
        )
        consensus.to_csv(OUT / "reverse_extension_consensus_top_pairs.csv", index=False)

    report_lines = [
        "# Reverse consistency extension report",
        "",
        f"Output directory: {OUT}",
        f"Analysis units: {inventory.shape[0]}",
        f"Repeat-level rows: {by_repeat.shape[0]}",
        "",
        "## Source groups",
    ]
    for source_group, sub in inventory.groupby("source_group"):
        report_lines.append(f"- {source_group}: {sub.shape[0]} dataset-condition-edge-set units")
    report_lines.extend(["", "## Metric overview"])
    report_lines.append("```csv")
    report_lines.append(overview.to_csv(index=False).strip())
    report_lines.append("```")
    (OUT / "reverse_extension_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print(f"Wrote {OUT}")
    print(f"analysis_units={inventory.shape[0]}")
    print(f"repeat_rows={by_repeat.shape[0]}")
    print(overview.to_string(index=False))


if __name__ == "__main__":
    analyze()
