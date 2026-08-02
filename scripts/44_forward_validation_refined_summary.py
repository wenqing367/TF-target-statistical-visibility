"""
Build unified tables for the forward known-edge visibility validation.

This script does not rerun matching or pair metrics. It only combines the
strict refined matching outputs from the local datasets and the two added
server datasets into a clean result source for the first project objective.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
LOCAL_REFINED = ROOT / "results" / "formal_abs_association_refined_matching"
SERVER_ADDED = ROOT / "results" / "server_added_refined_matching"
OUT = ROOT / "results" / "forward_validation_refined"

RESULT_SOURCES = [
    {
        "source": "local_refined",
        "result_dir": LOCAL_REFINED,
        "summary": LOCAL_REFINED / "refined_formal_abs_summary.csv",
        "by_repeat": LOCAL_REFINED / "refined_formal_abs_by_repeat.csv",
        "qc_summary": LOCAL_REFINED / "refined_matching_qc_summary.csv",
        "qc_by_repeat": LOCAL_REFINED / "refined_matching_qc_by_repeat.csv",
    },
    {
        "source": "server_nygc",
        "result_dir": SERVER_ADDED / "nygc_multimodal_pbmc_refined_matching",
        "summary": SERVER_ADDED / "nygc_multimodal_pbmc_refined_matching" / "refined_formal_abs_summary.csv",
        "by_repeat": SERVER_ADDED / "nygc_multimodal_pbmc_refined_matching" / "refined_formal_abs_by_repeat.csv",
        "qc_summary": SERVER_ADDED / "nygc_multimodal_pbmc_refined_matching" / "refined_matching_qc_summary.csv",
        "qc_by_repeat": SERVER_ADDED / "nygc_multimodal_pbmc_refined_matching" / "refined_matching_qc_by_repeat.csv",
    },
    {
        "source": "server_gse126030",
        "result_dir": SERVER_ADDED / "gse126030_refined_matching",
        "summary": SERVER_ADDED / "gse126030_refined_matching" / "refined_formal_abs_summary.csv",
        "by_repeat": SERVER_ADDED / "gse126030_refined_matching" / "refined_formal_abs_by_repeat.csv",
        "qc_summary": SERVER_ADDED / "gse126030_refined_matching" / "refined_matching_qc_summary.csv",
        "qc_by_repeat": SERVER_ADDED / "gse126030_refined_matching" / "refined_matching_qc_by_repeat.csv",
    },
]

EXPECTED_METRICS = {
    "pearson",
    "spearman",
    "mutual_information",
    "coexpression_probability",
    "codetection_odds_ratio",
}


def read_source(path: Path, source: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    df.insert(0, "result_source", source)
    return df


def check_required_columns(df: pd.DataFrame, required: set[str], name: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def add_dataset_role(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    role = []
    for _, row in out.iterrows():
        dataset = row["dataset"]
        condition = row["condition"]
        if dataset == "NYGC_multimodal_PBMC":
            role.append("normal broad PBMC validation")
        elif dataset == "GSE126030_T_cells":
            role.append(f"T-cell tissue validation ({condition})")
        elif dataset == "PBMC10k":
            role.append("original resting PBMC validation")
        elif dataset == "Kang_IFN_beta":
            role.append("original IFN-beta PBMC condition validation")
        elif dataset == "GSE178429_IFN_gamma":
            role.append("original IFN-gamma PBMC condition validation")
        else:
            role.append("forward validation")
    out["dataset_role"] = role
    return out


def build_completion_table(summary: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    units = (
        summary.groupby(["dataset", "condition", "edge_set"], dropna=False)
        .agg(
            n_metrics=("metric", "nunique"),
            n_summary_rows=("metric", "size"),
            min_repeats=("n_repeats", "min"),
            max_repeats=("n_repeats", "max"),
        )
        .reset_index()
    )
    qc_cols = [
        "dataset",
        "condition",
        "edge_set",
        "original_positive_edges",
        "matched_positive_edges_common",
        "matching_coverage",
        "per_repeat_positive_edges",
        "per_repeat_negative_min",
        "per_repeat_negative_max",
        "internal_duplicate_total",
        "internal_duplicate_max",
        "failed_matches_total_final",
        "relaxed_fraction_max",
        "positive_negative_overlap_unique",
        "negative_overlap_with_known_union_unique",
        "smd_target_mean_counts_max_abs",
        "smd_target_detection_rate_max_abs",
    ]
    present = [c for c in qc_cols if c in qc.columns]
    merged = units.merge(qc[present], on=["dataset", "condition", "edge_set"], how="left")
    merged["has_all_five_metrics"] = merged["n_metrics"].eq(len(EXPECTED_METRICS))
    merged["has_10_repeats"] = merged["min_repeats"].eq(10) & merged["max_repeats"].eq(10)
    merged["has_1to1_repeats"] = merged["per_repeat_positive_edges"].eq(merged["per_repeat_negative_min"]) & merged[
        "per_repeat_positive_edges"
    ].eq(merged["per_repeat_negative_max"])
    merged["qc_zero_duplicate_overlap_leakage"] = (
        merged["internal_duplicate_total"].fillna(0).eq(0)
        & merged["positive_negative_overlap_unique"].fillna(0).eq(0)
        & merged["negative_overlap_with_known_union_unique"].fillna(0).eq(0)
    )
    merged["forward_validation_ready"] = (
        merged["has_all_five_metrics"]
        & merged["has_10_repeats"]
        & merged["has_1to1_repeats"]
        & merged["qc_zero_duplicate_overlap_leakage"]
    )
    return merged.sort_values(["dataset", "condition", "edge_set"]).reset_index(drop=True)


def build_key_result_table(summary: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "result_source",
        "dataset",
        "condition",
        "dataset_role",
        "edge_set",
        "metric",
        "n_repeats",
        "n_positive_mean",
        "n_negative_mean",
        "random_auprc_baseline",
        "auroc_mean",
        "auroc_sd",
        "auroc_ci95_low",
        "auroc_ci95_high",
        "auprc_mean",
        "auprc_sd",
        "auprc_ci95_low",
        "auprc_ci95_high",
        "enrichment_at_5pct_mean",
        "median_difference_mean",
        "cliffs_delta_mean",
    ]
    present = [c for c in keep if c in summary.columns]
    out = summary[present].copy()
    out["auroc_delta_vs_0_5"] = out["auroc_mean"] - 0.5
    out["auprc_delta_vs_0_5"] = out["auprc_mean"] - 0.5
    return out.sort_values(["dataset", "condition", "edge_set", "metric"]).reset_index(drop=True)


def write_report(completion: pd.DataFrame, key_results: pd.DataFrame, qc: pd.DataFrame) -> None:
    n_units = len(completion)
    ready = int(completion["forward_validation_ready"].sum())
    n_rows = len(key_results)
    n_datasets = key_results[["dataset", "condition"]].drop_duplicates().shape[0]
    metric_counts = key_results.groupby("metric").size().to_dict()
    ps = key_results[key_results["metric"].isin(["pearson", "spearman"])]
    ps_positive = int((ps["auprc_mean"] > 0.5).sum())
    ps_total = len(ps)
    lines = [
        "# Forward Validation Refined Summary",
        "",
        "Scope: known TF-target edges versus strict expression/detection matched background edges.",
        "Pearson and Spearman are absolute correlations. No matching or metric recomputation is performed here.",
        "",
        f"- Analysis units: {n_units}",
        f"- Units passing completion/QC checks: {ready}/{n_units}",
        f"- Dataset-condition groups: {n_datasets}",
        f"- Metric result rows: {n_rows}",
        f"- Metric row counts: {metric_counts}",
        f"- Absolute Pearson/Spearman rows with AUPRC > 0.5: {ps_positive}/{ps_total}",
        "",
        "Important dataset note: GSE126030 is represented by tissue-specific T-cell h5ad files. The server files do not contain condition/stim/sample metadata, so these results are tissue validation results, not anti-CD3/CD28 stimulation comparisons.",
        "",
        "Output tables:",
        "- forward_refined_key_results.csv",
        "- forward_refined_qc_summary.csv",
        "- forward_refined_by_repeat.csv",
        "- forward_refined_qc_by_repeat.csv",
        "- forward_refined_completion_check.csv",
    ]
    (OUT / "forward_refined_summary_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = pd.concat([read_source(s["summary"], s["source"]) for s in RESULT_SOURCES], ignore_index=True)
    by_repeat = pd.concat([read_source(s["by_repeat"], s["source"]) for s in RESULT_SOURCES], ignore_index=True)
    qc = pd.concat([read_source(s["qc_summary"], s["source"]) for s in RESULT_SOURCES], ignore_index=True)
    qc_by_repeat = pd.concat([read_source(s["qc_by_repeat"], s["source"]) for s in RESULT_SOURCES], ignore_index=True)

    check_required_columns(
        summary,
        {"dataset", "condition", "edge_set", "metric", "n_repeats", "auroc_mean", "auprc_mean"},
        "summary",
    )
    check_required_columns(
        qc,
        {
            "dataset",
            "condition",
            "edge_set",
            "matching_coverage",
            "positive_negative_overlap_unique",
            "negative_overlap_with_known_union_unique",
        },
        "qc summary",
    )

    observed_metrics = set(summary["metric"].unique())
    if observed_metrics != EXPECTED_METRICS:
        raise ValueError(f"Unexpected metric set: {sorted(observed_metrics)}")

    summary = add_dataset_role(summary)
    key_results = build_key_result_table(summary)
    completion = build_completion_table(summary, qc)

    key_results.to_csv(OUT / "forward_refined_key_results.csv", index=False)
    summary.to_csv(OUT / "forward_refined_full_summary.csv", index=False)
    by_repeat.to_csv(OUT / "forward_refined_by_repeat.csv", index=False)
    qc.to_csv(OUT / "forward_refined_qc_summary.csv", index=False)
    qc_by_repeat.to_csv(OUT / "forward_refined_qc_by_repeat.csv", index=False)
    completion.to_csv(OUT / "forward_refined_completion_check.csv", index=False)
    write_report(completion, key_results, qc)
    print(f"[write] {OUT}")
    print(f"[ready] {int(completion['forward_validation_ready'].sum())}/{len(completion)} units")


if __name__ == "__main__":
    main()
