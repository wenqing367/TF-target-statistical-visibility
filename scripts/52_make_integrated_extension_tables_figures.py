# -*- coding: utf-8 -*-
"""Build integrated paper tables and figures for original plus extension analyses.

Outputs are written to:
  results/integrated_paper_outputs/
  figures/integrated/

The script does not recompute matched background pairs or pair-level metrics.
It only summarizes included refined-result tables and, when available, builds
Hallmark-wide functional-coherence summaries from stored top-pair results.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[1]
RESULT_OUT = ROOT / "results" / "integrated_paper_outputs"
FIG_OUT = ROOT / "figures" / "integrated"

METRIC_ORDER = [
    "pearson",
    "spearman",
    "mutual_information",
    "coexpression_probability",
    "codetection_odds_ratio",
]
METRIC_LABELS = {
    "pearson": "|Pearson|",
    "spearman": "|Spearman|",
    "mutual_information": "MI",
    "coexpression_probability": "Co-expression",
    "codetection_odds_ratio": "Co-detection OR",
}
PBMC_DATASET_LABELS = {
    "PBMC10k": "PBMC10k",
    "Kang_IFN_beta": "Kang\nIFN-beta\nPBMC",
    "GSE178429_IFN_gamma": "GSE178429\nIFN-gamma\nPBMC",
    "NYGC_multimodal_PBMC": "NYGC\nPBMC",
}
ADULT_TISSUE_LABELS = {
    "AdultColon": "Adult\ncolon",
    "AdultKidney": "Adult\nkidney",
    "AdultLiver": "Adult\nliver",
    "AdultLung": "Adult\nlung",
}
TISSUE_CONTEXT_ORDER = [
    "PBMC10k",
    "Kang\nIFN-beta\nPBMC",
    "GSE178429\nIFN-gamma\nPBMC",
    "NYGC\nPBMC",
    "Adult\ncolon",
    "Adult\nkidney",
    "Adult\nliver",
    "Adult\nlung",
]
CELLTYPE_CONTEXT_ORDER = [
    "AT2 cell",
    "Endothelial cell",
    "Enterocyte",
    "Epithelial cell",
    "Fibroblast",
    "Loop of Henle",
    "Smooth muscle cell",
]
REVISED_CONTEXT_LEVEL_ORDER = ["Tissue contexts", "Adult cell types"]
REVISED_CONTEXT_LABELS = {
    "Tissue contexts": "PBMC + adult\ntissue contexts",
    "Adult cell types": "Adult\ncell types",
}
HALLMARK_GROUP_ORDER = ["PBMC datasets", "Adult tissues", "Adult cell types"]
HALLMARK_GROUP_LABELS = {
    "PBMC datasets": "A. PBMC datasets",
    "Adult tissues": "B. Adult tissues",
    "Adult cell types": "C. Adult cell types",
}
HALLMARK_GROUP_COLORS = {
    "PBMC datasets": "#6B8FB3",
    "Adult tissues": "#D18F45",
    "Adult cell types": "#8AAE92",
}
HALLMARK_DISPLAY_LABELS = {
    "Allograft Rejection": "Allograft Rejection",
    "Apoptosis": "Apoptosis",
    "Complement": "Complement",
    "Epithelial Mesenchymal Transition": "Epithelial-Mesenchymal Transition",
    "Hypoxia": "Hypoxia",
    "Il2 Stat5 Signaling": "IL-2/STAT5 Signaling",
    "Il6 Jak Stat3 Signaling": "IL-6/JAK/STAT3 Signaling",
    "Inflammatory Response": "Inflammatory Response",
    "Interferon Alpha Response": "IFN-α response",
    "Interferon Gamma Response": "IFN-γ response",
    "Kras Signaling Up": "KRAS Signaling Up",
    "Mtorc1 Signaling": "mTORC1 Signaling",
    "P53 Pathway": "p53 Pathway",
    "Tnfa Signaling Via Nfkb": "TNFα signaling via NF-κB",
}
COLORS = {
    "pearson": "#5F83A6",
    "spearman": "#7FA98F",
    "mutual_information": "#AEB8C2",
    "coexpression_probability": "#C8B783",
    "codetection_odds_ratio": "#A184A5",
}


def setup() -> None:
    RESULT_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 8,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
        }
    )


def ci95(values: pd.Series) -> float:
    values = values.dropna()
    if len(values) <= 1:
        return 0.0
    return float(1.96 * values.std(ddof=1) / np.sqrt(len(values)))


def sd_sample(values: pd.Series) -> float:
    values = values.dropna()
    if len(values) <= 1:
        return 0.0
    return float(values.std(ddof=1))


def save_figure(fig: plt.Figure, name: str) -> None:
    png_path = FIG_OUT / f"{name}.png"
    svg_path = FIG_OUT / f"{name}.svg"
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    print(png_path)
    print(svg_path)
    plt.close(fig)


def standardize_forward(df: pd.DataFrame, group: str) -> pd.DataFrame:
    out = df.copy()
    out["analysis_group"] = group
    return out


def load_forward_tables() -> pd.DataFrame:
    frames = []
    frames.append(
        standardize_forward(
            pd.read_csv(ROOT / "results" / "forward_validation_refined" / "forward_refined_full_summary.csv"),
            "Original dataset contexts",
        )
    )
    frames.append(
        standardize_forward(
            pd.read_csv(
                ROOT
                / "results"
                / "extension_forward_visibility"
                / "hcl_adult_tissue"
                / "zju_hcl_adult_tissue_refined_formal_abs_summary_combined.csv"
            ),
            "Adult tissues",
        )
    )
    frames.append(
        standardize_forward(
            pd.read_csv(
                ROOT
                / "results"
                / "extension_forward_visibility"
                / "hcl_adult_celltype"
                / "zju_hcl_adult_celltype_refined_formal_abs_summary_combined.csv"
            ),
            "Adult cell types",
        )
    )
    return pd.concat(frames, ignore_index=True, sort=False)


def summarize_forward(forward: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_summary = (
        forward.groupby(["analysis_group", "metric"], dropna=False)
        .agg(
            analysis_units=("auprc_mean", "size"),
            mean_auprc=("auprc_mean", "mean"),
            ci95_auprc=("auprc_mean", ci95),
            mean_auroc=("auroc_mean", "mean"),
            ci95_auroc=("auroc_mean", ci95),
            units_auprc_above_0_5=("auprc_mean", lambda x: int((x > 0.5).sum())),
            mean_cliffs_delta=("cliffs_delta_mean", "mean"),
            median_matched_positives=("n_positive_mean", "median"),
        )
        .reset_index()
    )
    context_summary = (
        forward.groupby(["analysis_group", "dataset", "condition", "edge_set"], dropna=False)
        .agg(
            metrics=("metric", "nunique"),
            median_matched_positives=("n_positive_mean", "median"),
            mean_auprc=("auprc_mean", "mean"),
            mean_auroc=("auroc_mean", "mean"),
            mean_cliffs_delta=("cliffs_delta_mean", "mean"),
        )
        .reset_index()
    )
    metric_summary.to_csv(RESULT_OUT / "integrated_forward_metric_summary.csv", index=False)
    context_summary.to_csv(RESULT_OUT / "integrated_forward_context_summary.csv", index=False)
    return metric_summary, context_summary


def load_reverse_tables() -> pd.DataFrame:
    original = pd.read_csv(
        ROOT / "results" / "reverse_validation_from_refined_pairs" / "reverse_refined_known_overlap_summary.csv"
    )
    original = original.rename(
        columns={
            "n_repeats": "repeats",
            "known_positive_fraction_mean": "mean_curated_positive_fraction",
            "known_positive_fraction_sd": "sd_curated_positive_fraction",
            "delta_vs_balanced_baseline_mean": "mean_delta_vs_0_5",
            "top_n_mean": "mean_actual_top_n",
        }
    )
    original["analysis_group"] = "Original dataset contexts"
    extension = pd.read_csv(
        ROOT / "results" / "extension_reverse_consistency" / "reverse_extension_known_fraction_summary.csv"
    )
    extension = extension[extension["source_group"].isin(["hcl_adult_tissue", "hcl_adult_celltype"])].copy()
    group_map = {
        "hcl_adult_tissue": "Adult tissues",
        "hcl_adult_celltype": "Adult cell types",
    }
    extension["analysis_group"] = extension["source_group"].map(group_map)
    return pd.concat([original, extension], ignore_index=True, sort=False)


def summarize_reverse(reverse: pd.DataFrame) -> pd.DataFrame:
    summary = (
        reverse.groupby(["analysis_group", "metric", "top_level"], dropna=False)
        .agg(
            analysis_units=("mean_curated_positive_fraction", "size"),
            mean_curated_positive_fraction=("mean_curated_positive_fraction", "mean"),
            ci95_curated_positive_fraction=("mean_curated_positive_fraction", ci95),
            units_mean_above_0_5=("mean_curated_positive_fraction", lambda x: int((x > 0.5).sum())),
            median_actual_top_n=("mean_actual_top_n", "median"),
            mean_delta_vs_0_5=("mean_delta_vs_0_5", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(RESULT_OUT / "integrated_reverse_metric_summary.csv", index=False)
    reverse.to_csv(RESULT_OUT / "integrated_reverse_unit_summary.csv", index=False)
    return summary


def write_supplementary_workbook(
    forward_metric_summary: pd.DataFrame,
    forward_context_summary: pd.DataFrame,
    reverse_metric_summary: pd.DataFrame,
    reverse_unit_summary: pd.DataFrame,
    hallmark_overview: pd.DataFrame,
) -> None:
    """Write the compact manuscript Supplementary Tables S1-S5."""

    qc_paths = [
        ROOT / "results" / "forward_validation_refined" / "forward_refined_qc_summary.csv",
        ROOT
        / "results"
        / "extension_forward_visibility"
        / "hcl_adult_tissue"
        / "zju_hcl_adult_tissue_refined_matching_qc_summary_combined.csv",
        ROOT
        / "results"
        / "extension_forward_visibility"
        / "hcl_adult_celltype"
        / "zju_hcl_adult_celltype_refined_matching_qc_summary_combined.csv",
    ]
    qc = pd.concat([pd.read_csv(path) for path in qc_paths], ignore_index=True, sort=False)
    qc = add_revised_context_columns(qc)

    edge_labels = {
        "trrust": "TRRUST",
        "dorothea_ab": "DoRothEA A/B",
        "trrust_dorothea_intersection": "intersection",
        "trrust_dorothea_union": "union",
    }
    context_order = TISSUE_CONTEXT_ORDER + CELLTYPE_CONTEXT_ORDER
    context_rank = {label: i for i, label in enumerate(context_order)}

    coverage_rows = []
    for (context_level, context_label), sub in qc.groupby(["context_level", "context_label"], sort=False):
        dataset = str(sub["dataset"].iloc[0])
        if dataset in PBMC_DATASET_LABELS:
            biological_group = "PBMC datasets"
            conditions = "; ".join(dict.fromkeys(sub["condition"].astype(str)))
        elif dataset == "ZJU_HCL_adult_tissue":
            biological_group = "Adult tissues"
            conditions = "all"
        else:
            biological_group = "Adult cell types"
            conditions = "all"
        edge_sets = [edge_labels[name] for name in edge_labels if name in set(sub["edge_set"])]
        retained_min = int(sub["matched_positive_edges_common"].min())
        retained_max = int(sub["matched_positive_edges_common"].max())
        coverage_rows.append(
            {
                "Biological group": biological_group,
                "Dataset / final context": str(context_label).replace("\n", " "),
                "Condition(s)": conditions,
                "Available edge sets": "; ".join(edge_sets),
                "Analysis units": int(len(sub)),
                "Retained positives": f"{retained_min}-{retained_max}",
                "Matching coverage": (
                    f"{100 * sub['matching_coverage'].min():.1f}%-"
                    f"{100 * sub['matching_coverage'].max():.1f}%"
                ),
                "Repeats": int(sub["repeats"].max()),
                "_order": context_rank[str(context_label)],
            }
        )
    table_s1 = pd.DataFrame(coverage_rows).sort_values("_order").drop(columns="_order")

    table_s2 = pd.DataFrame(
        [
            ("Condition × edge-set analysis units", len(qc)),
            ("Matched-background repeats per unit", int(qc["repeats"].max())),
            (
                "Retained positive edges per unit",
                f"{int(qc['matched_positive_edges_common'].min())}-{int(qc['matched_positive_edges_common'].max())}",
            ),
            (
                "Matching coverage across units",
                f"{100 * qc['matching_coverage'].min():.1f}%-{100 * qc['matching_coverage'].max():.1f}%",
            ),
            (
                "Equal positive/background counts in every repeat",
                f"{int((qc['per_repeat_positive_edges'].eq(qc['per_repeat_negative_min']) & qc['per_repeat_positive_edges'].eq(qc['per_repeat_negative_max'])).sum())}/{len(qc)}",
            ),
            ("Duplicated background pairs", int(qc["internal_duplicate_total"].max())),
            ("Positive-background overlap", int(qc["positive_negative_overlap_unique"].max())),
            (
                "Background leakage into curated-edge union",
                int(qc["negative_overlap_with_known_union_unique"].max()),
            ),
            ("Maximum absolute SMD, target mean expression", f"{qc['smd_target_mean_counts_max_abs'].max():.3f}"),
            ("Maximum absolute SMD, target detection rate", f"{qc['smd_target_detection_rate_max_abs'].max():.3f}"),
        ],
        columns=["QC item", "Value"],
    )

    context_group_labels = {
        "Tissue contexts": "PBMC + adult tissues",
        "Adult cell types": "Adult cell types",
    }

    table_s3 = forward_metric_summary.copy()
    table_s3["Context group"] = table_s3["context_level"].map(context_group_labels)
    table_s3["Metric"] = table_s3["metric"].map(METRIC_LABELS)
    table_s3["_group"] = table_s3["context_level"].map({"Tissue contexts": 0, "Adult cell types": 1})
    table_s3["_metric"] = table_s3["metric"].map({name: i for i, name in enumerate(METRIC_ORDER)})
    table_s3 = table_s3.sort_values(["_group", "_metric"])
    table_s3 = table_s3.assign(
        **{
            "Final contexts": table_s3["analysis_units"].astype(int),
            "Mean AUPRC": table_s3["mean_auprc"].round(3),
            "Mean AUROC": table_s3["mean_auroc"].round(3),
            "Mean Cliff's delta": table_s3["mean_cliffs_delta"].round(3),
            "AUPRC > 0.5 contexts": (
                table_s3["units_auprc_above_0_5"].astype(int).astype(str)
                + "/"
                + table_s3["analysis_units"].astype(int).astype(str)
            ),
        }
    )[
        [
            "Context group",
            "Metric",
            "Final contexts",
            "Mean AUPRC",
            "Mean AUROC",
            "Mean Cliff's delta",
            "AUPRC > 0.5 contexts",
        ]
    ]

    reverse = reverse_metric_summary.copy()
    reverse["Context group"] = reverse["context_level"].map(context_group_labels)
    reverse["Metric"] = reverse["metric"].map(METRIC_LABELS)
    reverse["_group"] = reverse["context_level"].map({"Tissue contexts": 0, "Adult cell types": 1})
    reverse["_metric"] = reverse["metric"].map({name: i for i, name in enumerate(METRIC_ORDER)})

    table_s4a = reverse[reverse["top_level"].eq("top_100")].sort_values(["_group", "_metric"])
    table_s4a = table_s4a.assign(
        **{
            "Valid final contexts": table_s4a["analysis_units"].astype(int),
            "Mean Top100 fraction": table_s4a["mean_curated_positive_fraction"].round(3),
            "Contexts > 0.5": (
                table_s4a["units_mean_above_0_5"].astype(int).astype(str)
                + "/"
                + table_s4a["analysis_units"].astype(int).astype(str)
            ),
            "Delta vs 0.5": table_s4a["mean_delta_vs_0_5"].round(3),
        }
    )[
        ["Context group", "Metric", "Valid final contexts", "Mean Top100 fraction", "Contexts > 0.5", "Delta vs 0.5"]
    ]

    table_s4b = reverse[reverse["top_level"].isin(["top_500", "top_5pct"])].copy()
    table_s4b["_sensitivity"] = table_s4b["top_level"].map({"top_500": 0, "top_5pct": 1})
    table_s4b["Sensitivity"] = table_s4b["top_level"].map({"top_500": "Strict Top500", "top_5pct": "top-5%"})
    table_s4b = table_s4b.sort_values(["_sensitivity", "_group", "_metric"])
    table_s4b = table_s4b.assign(
        **{
            "Valid final contexts": table_s4b["analysis_units"].astype(int),
            "Mean curated-positive fraction": table_s4b["mean_curated_positive_fraction"].round(3),
            "Contexts > 0.5": (
                table_s4b["units_mean_above_0_5"].astype(int).astype(str)
                + "/"
                + table_s4b["analysis_units"].astype(int).astype(str)
            ),
            "Delta vs 0.5": table_s4b["mean_delta_vs_0_5"].round(3),
        }
    )[
        [
            "Sensitivity",
            "Context group",
            "Metric",
            "Valid final contexts",
            "Mean curated-positive fraction",
            "Contexts > 0.5",
            "Delta vs 0.5",
        ]
    ]

    hallmark_rows = []
    for group in HALLMARK_GROUP_ORDER:
        sub = hallmark_overview[
            hallmark_overview["context_group"].eq(group) & hallmark_overview["top_level"].eq("top_100")
        ].copy()
        sub = sub.sort_values(["support_fraction", "median_fold_enrichment"], ascending=False).head(7)
        for _, row in sub.iterrows():
            hallmark_rows.append(
                {
                    "Group": group,
                    "Pathway": HALLMARK_DISPLAY_LABELS.get(row["pathway"], row["pathway"]),
                    "Median FE": round(float(row["median_fold_enrichment"]), 3),
                    "Supported combinations": int(row["combinations_fdr_in_at_least_5_repeats"]),
                    "Analyzable combinations": int(row["combinations"]),
                    "Support fraction": round(float(row["support_fraction"]), 3),
                }
            )
    table_s5 = pd.DataFrame(hallmark_rows)

    sheets = {
        "S1_coverage": table_s1,
        "S2_matching_QC": table_s2,
        "S3_forward": table_s3,
        "S4A_fixed_Top100": table_s4a,
        "S4B_sensitivity": table_s4b,
        "S5_Hallmark_Figure4": table_s5,
    }
    out_path = RESULT_OUT / "integrated_supplementary_tables.xlsx"
    try:
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            for sheet_name, table in sheets.items():
                table.to_excel(writer, sheet_name=sheet_name, index=False)
            for sheet in writer.book.worksheets:
                sheet.freeze_panes = "A2"
                for col_cells in sheet.columns:
                    max_len = max(len(str(cell.value)) if cell.value is not None else 0 for cell in col_cells)
                    sheet.column_dimensions[col_cells[0].column_letter].width = min(max(max_len + 2, 10), 34)
    except ImportError:
        print("[warn] openpyxl is unavailable; skipped integrated_supplementary_tables.xlsx")
        return
    print(f"[write] {out_path}")


def standardize_gene_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper()


def standardize_edges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out["target"] = out["target"].astype(str).str.upper()
    return out


def prefix(dataset: str, condition: str, edge_set: str | None = None) -> str:
    parts = [dataset, condition]
    if edge_set is not None:
        parts.append(edge_set)
    return "_".join(parts).replace("/", "_")


def load_hallmark_sets() -> pd.DataFrame:
    path = ROOT / "results" / "extension_hallmark_coherence" / "hallmark_full_gene_sets_used.csv"
    df = pd.read_csv(path)
    if "pathway" not in df.columns or "gene" not in df.columns:
        raise ValueError(f"Invalid Hallmark file: {path}")
    df["gene"] = standardize_gene_series(df["gene"])
    return df[["collection", "pathway", "gene"]].drop_duplicates()


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


def source_dir_for(result_source: str) -> Path:
    source_dirs = {
        "local_refined": ROOT / "results" / "formal_abs_association_refined_matching",
        "server_nygc": ROOT / "results" / "server_added_refined_matching" / "nygc_multimodal_pbmc_refined_matching",
    }
    return source_dirs[result_source]


def build_original_universe_lookup(units: pd.DataFrame) -> dict[tuple, set[str]]:
    lookup: dict[tuple, set[str]] = {}
    for row in units.itertuples(index=False):
        result_source = row.result_source
        dataset = row.dataset
        condition = row.condition
        edge_set = row.edge_set
        source_dir = source_dir_for(result_source)
        pfx = prefix(dataset, condition, edge_set)
        pos_path = source_dir / f"{pfx}_refined_positives.csv"
        neg_path = source_dir / f"{pfx}_refined_negatives.csv"
        if not pos_path.exists() or not neg_path.exists():
            continue
        positives = pd.read_csv(pos_path, usecols=["target"])
        negatives = pd.read_csv(neg_path, usecols=["repeat", "target"])
        pos_targets = set(standardize_gene_series(positives["target"]))
        negatives["target"] = standardize_gene_series(negatives["target"])
        for repeat, group in negatives.groupby("repeat"):
            lookup[(result_source, dataset, condition, edge_set, int(repeat))] = pos_targets | set(group["target"])
    return lookup


def run_hallmark_tests_for_original() -> pd.DataFrame:
    out_path = RESULT_OUT / "original_hallmark_full_summary.csv"
    overview_path = RESULT_OUT / "original_hallmark_full_overview.csv"
    if out_path.exists() and overview_path.exists():
        return pd.read_csv(overview_path)

    top_path = ROOT / "results" / "reverse_validation_from_refined_pairs" / "reverse_refined_top_pairs_by_repeat.csv"
    if not top_path.exists():
        print("[skip] original Hallmark full analysis requires reverse_refined_top_pairs_by_repeat.csv")
        return pd.DataFrame()

    top_pairs = pd.read_csv(top_path)
    top_pairs = top_pairs[top_pairs["top_level"].isin(["top_100", "top_500", "top_5pct"])].copy()
    top_pairs = top_pairs[top_pairs["dataset"].isin(PBMC_DATASET_LABELS)].copy()
    units = top_pairs[["result_source", "dataset", "condition", "edge_set"]].drop_duplicates()
    universe_lookup = build_original_universe_lookup(units)
    pathways = load_hallmark_sets()
    pathway_sets = {
        key: set(group["gene"])
        for key, group in pathways.groupby(["collection", "pathway"], sort=False)
    }
    rows = []
    group_cols = ["result_source", "dataset", "condition", "edge_set", "repeat", "metric", "top_level"]
    for key, group in top_pairs.groupby(group_cols, sort=False):
        result_source, dataset, condition, edge_set, repeat, metric, top_level = key
        universe = universe_lookup.get((result_source, dataset, condition, edge_set, int(repeat)))
        if not universe:
            continue
        selected = group[group["label"].eq(1)]
        selected_targets = set(standardize_gene_series(selected["target"]))
        if not selected_targets:
            continue
        m = len(universe)
        n = len(selected_targets)
        family_rows = []
        for (collection, pathway), genes in pathway_sets.items():
            geneset = genes & universe
            if not geneset:
                continue
            overlap = selected_targets & geneset
            k = len(overlap)
            fold = (k / n) / (len(geneset) / m) if n and geneset else np.nan
            pvalue = hypergeom.sf(k - 1, m, len(geneset), n) if k > 0 else 1.0
            family_rows.append(
                {
                    "analysis_group": "Original dataset contexts",
                    "result_source": result_source,
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
                    "fold_enrichment": fold,
                    "pvalue": pvalue,
                }
            )
        if family_rows:
            fam = pd.DataFrame(family_rows)
            fam["bh_fdr"] = bh_fdr(fam["pvalue"])
            rows.append(fam)
    tests = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if tests.empty:
        return pd.DataFrame()
    summary = (
        tests.groupby(
            ["analysis_group", "dataset", "condition", "edge_set", "metric", "top_level", "collection", "pathway"],
            dropna=False,
        )
        .agg(
            repeats_tested=("repeat", "nunique"),
            repeats_with_overlap=("overlap_genes", lambda x: int((x > 0).sum())),
            repeats_fdr_lt_0_05=("bh_fdr", lambda x: int((x < 0.05).sum())),
            median_selected_targets=("selected_target_count", "median"),
            median_overlap_genes=("overlap_genes", "median"),
            median_fold_enrichment=("fold_enrichment", "median"),
            best_bh_fdr=("bh_fdr", "min"),
        )
        .reset_index()
    )
    overview = (
        summary.groupby(["analysis_group", "pathway", "top_level"], dropna=False)
        .agg(
            combinations=("metric", "count"),
            combinations_fdr_in_at_least_5_repeats=("repeats_fdr_lt_0_05", lambda x: int((x >= 5).sum())),
            median_fold_enrichment=("median_fold_enrichment", "median"),
            best_bh_fdr=("best_bh_fdr", "min"),
        )
        .reset_index()
    )
    summary.to_csv(out_path, index=False)
    overview.to_csv(overview_path, index=False)
    return overview


def load_hallmark_overview() -> pd.DataFrame:
    original = run_hallmark_tests_for_original()
    extension = pd.read_csv(ROOT / "results" / "extension_hallmark_coherence" / "hallmark_extension_overview.csv")
    extension = extension[extension["source_group"].isin(["hcl_adult_tissue", "hcl_adult_celltype"])].copy()
    group_map = {
        "hcl_adult_tissue": "Adult tissues",
        "hcl_adult_celltype": "Adult cell types",
    }
    extension["analysis_group"] = extension["source_group"].map(group_map)
    keep = [
        "analysis_group",
        "pathway",
        "top_level",
        "combinations",
        "combinations_fdr_in_at_least_5_repeats",
        "median_fold_enrichment",
        "best_bh_fdr",
    ]
    frames = []
    if not original.empty:
        frames.append(original[keep])
    frames.append(extension[keep])
    out = pd.concat(frames, ignore_index=True)
    out["support_fraction"] = out["combinations_fdr_in_at_least_5_repeats"] / out["combinations"]
    out.to_csv(RESULT_OUT / "integrated_hallmark_overview.csv", index=False)
    return out


def add_revised_context_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the revised manuscript scope and add display-level context labels."""
    out = df.copy()
    out["context_level"] = pd.NA
    out["context_label"] = pd.NA

    pbmc_mask = out["dataset"].isin(PBMC_DATASET_LABELS)
    out.loc[pbmc_mask, "context_level"] = "Tissue contexts"
    out.loc[pbmc_mask, "context_label"] = out.loc[pbmc_mask, "dataset"].map(PBMC_DATASET_LABELS)

    adult_tissue_mask = out["dataset"].eq("ZJU_HCL_adult_tissue")
    out.loc[adult_tissue_mask, "context_level"] = "Tissue contexts"
    out.loc[adult_tissue_mask, "context_label"] = out.loc[adult_tissue_mask, "condition"].map(ADULT_TISSUE_LABELS)

    adult_cell_mask = out["dataset"].eq("ZJU_HCL_adult_celltype")
    out.loc[adult_cell_mask, "context_level"] = "Adult cell types"
    out.loc[adult_cell_mask, "context_label"] = out.loc[adult_cell_mask, "condition"]

    keep = out["context_level"].notna() & out["context_label"].notna()
    return out.loc[keep].copy()


def summarize_forward_revised(forward: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    revised = add_revised_context_columns(forward)
    context_summary = (
        revised.groupby(["context_level", "context_label", "metric"], dropna=False)
        .agg(
            analysis_units=("auprc_mean", "size"),
            mean_auprc=("auprc_mean", "mean"),
            ci95_auprc=("auprc_mean", ci95),
            mean_auroc=("auroc_mean", "mean"),
            ci95_auroc=("auroc_mean", ci95),
            units_auprc_above_0_5=("auprc_mean", lambda x: int((x > 0.5).sum())),
            mean_cliffs_delta=("cliffs_delta_mean", "mean"),
            median_matched_positives=("n_positive_mean", "median"),
        )
        .reset_index()
    )
    metric_summary = (
        context_summary.groupby(["context_level", "metric"], dropna=False)
        .agg(
            analysis_units=("mean_auprc", "size"),
            mean_auprc=("mean_auprc", "mean"),
            sd_auprc=("mean_auprc", sd_sample),
            ci95_auprc=("mean_auprc", ci95),
            mean_auroc=("mean_auroc", "mean"),
            sd_auroc=("mean_auroc", sd_sample),
            ci95_auroc=("mean_auroc", ci95),
            units_auprc_above_0_5=("mean_auprc", lambda x: int((x > 0.5).sum())),
            mean_cliffs_delta=("mean_cliffs_delta", "mean"),
            sd_cliffs_delta=("mean_cliffs_delta", sd_sample),
            median_matched_positives=("median_matched_positives", "median"),
        )
        .reset_index()
    )
    context_summary.to_csv(RESULT_OUT / "revised_forward_context_summary.csv", index=False)
    metric_summary.to_csv(RESULT_OUT / "revised_forward_metric_summary.csv", index=False)
    revised.to_csv(RESULT_OUT / "revised_forward_analysis_units.csv", index=False)
    context_summary.to_csv(RESULT_OUT / "integrated_forward_context_summary.csv", index=False)
    metric_summary.to_csv(RESULT_OUT / "integrated_forward_metric_summary.csv", index=False)
    return context_summary, metric_summary


def summarize_reverse_revised(reverse: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    revised = add_revised_context_columns(reverse)

    def incomplete_fixed_threshold(top_level: str, threshold: int) -> pd.Series:
        at_level = revised["top_level"].eq(top_level)
        incomplete = at_level & revised["mean_actual_top_n"].lt(threshold)
        if "mean_total_pairs_ranked" in revised.columns:
            incomplete |= (
                at_level
                & revised["mean_total_pairs_ranked"].notna()
                & revised["mean_total_pairs_ranked"].lt(threshold)
            )
        if "complete_threshold_repeats" in revised.columns and "repeats" in revised.columns:
            incomplete |= (
                at_level
                & revised["complete_threshold_repeats"].notna()
                & revised["repeats"].notna()
                & revised["complete_threshold_repeats"].lt(revised["repeats"])
            )
        return incomplete

    incomplete_top100 = incomplete_fixed_threshold("top_100", 100)
    incomplete_top500 = incomplete_fixed_threshold("top_500", 500)
    revised["fixed_top100_valid"] = ~incomplete_top100
    revised["strict_top500_valid"] = ~incomplete_top500
    invalid_fixed_threshold = incomplete_top100 | incomplete_top500
    for col in [
        "mean_curated_positive_fraction",
        "mean_delta_vs_0_5",
        "repeats_above_0_5",
        "ci95_curated_positive_fraction",
        "min_curated_positive_fraction",
        "max_curated_positive_fraction",
    ]:
        if col in revised.columns:
            revised.loc[invalid_fixed_threshold, col] = np.nan

    context_summary = (
        revised.groupby(["context_level", "context_label", "metric", "top_level"], dropna=False)
        .agg(
            analysis_units=("mean_curated_positive_fraction", "count"),
            mean_curated_positive_fraction=("mean_curated_positive_fraction", "mean"),
            ci95_curated_positive_fraction=("mean_curated_positive_fraction", ci95),
            units_mean_above_0_5=("mean_curated_positive_fraction", lambda x: int((x > 0.5).sum())),
            median_actual_top_n=("mean_actual_top_n", "median"),
            mean_delta_vs_0_5=("mean_delta_vs_0_5", "mean"),
        )
        .reset_index()
    )
    metric_summary = (
        context_summary.groupby(["context_level", "metric", "top_level"], dropna=False)
        .agg(
            analysis_units=("mean_curated_positive_fraction", "count"),
            mean_curated_positive_fraction=("mean_curated_positive_fraction", "mean"),
            sd_curated_positive_fraction=("mean_curated_positive_fraction", sd_sample),
            ci95_curated_positive_fraction=("mean_curated_positive_fraction", ci95),
            units_mean_above_0_5=("mean_curated_positive_fraction", lambda x: int((x > 0.5).sum())),
            median_actual_top_n=("median_actual_top_n", "median"),
            mean_delta_vs_0_5=("mean_delta_vs_0_5", "mean"),
        )
        .reset_index()
    )
    context_summary.to_csv(RESULT_OUT / "revised_reverse_context_summary.csv", index=False)
    metric_summary.to_csv(RESULT_OUT / "revised_reverse_metric_summary.csv", index=False)
    revised.to_csv(RESULT_OUT / "revised_reverse_analysis_units.csv", index=False)
    context_summary.to_csv(RESULT_OUT / "integrated_reverse_unit_summary.csv", index=False)
    metric_summary.to_csv(RESULT_OUT / "integrated_reverse_metric_summary.csv", index=False)
    return context_summary, metric_summary


def load_revised_hallmark_detail() -> pd.DataFrame:
    run_hallmark_tests_for_original()
    original = pd.read_csv(RESULT_OUT / "original_hallmark_full_summary.csv")
    original["source_group"] = "original"
    extension = pd.read_csv(ROOT / "results" / "extension_hallmark_coherence" / "hallmark_extension_summary.csv")
    extension = extension[extension["source_group"].isin(["hcl_adult_tissue", "hcl_adult_celltype"])].copy()
    detail = pd.concat([original, extension], ignore_index=True, sort=False)
    detail = add_revised_context_columns(detail)
    return detail


def summarize_hallmark_revised() -> pd.DataFrame:
    detail = load_revised_hallmark_detail()
    detail["repeat_supported"] = detail["repeats_fdr_lt_0_05"] >= 5
    detail["context_group"] = pd.NA
    detail.loc[detail["dataset"].isin(PBMC_DATASET_LABELS), "context_group"] = "PBMC datasets"
    detail.loc[detail["dataset"].eq("ZJU_HCL_adult_tissue"), "context_group"] = "Adult tissues"
    detail.loc[detail["dataset"].eq("ZJU_HCL_adult_celltype"), "context_group"] = "Adult cell types"
    detail = detail[detail["context_group"].notna()].copy()
    overview = (
        detail.groupby(["context_group", "pathway", "top_level"], dropna=False)
        .agg(
            combinations=("repeat_supported", "size"),
            combinations_fdr_in_at_least_5_repeats=("repeat_supported", "sum"),
            median_fold_enrichment=("median_fold_enrichment", "median"),
            best_bh_fdr=("best_bh_fdr", "min"),
        )
        .reset_index()
    )
    overview["support_fraction"] = overview["combinations_fdr_in_at_least_5_repeats"] / overview["combinations"]
    overview["context_level"] = overview["context_group"]
    overview.to_csv(RESULT_OUT / "revised_hallmark_context_summary.csv", index=False)
    detail.to_csv(RESULT_OUT / "revised_hallmark_analysis_units.csv", index=False)
    overview.to_csv(RESULT_OUT / "integrated_hallmark_overview.csv", index=False)
    return overview


def add_value_labels(ax, bars, vals, errs=None, fmt="{:.3f}", dy=0.002) -> None:
    if errs is None:
        errs = np.zeros(len(vals))
    for bar, val, err in zip(bars, vals, errs):
        if pd.isna(val):
            continue
        va = "bottom" if val >= 0.5 else "top"
        y = val + float(err) + dy if val >= 0.5 else val - float(err) - dy
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y,
            fmt.format(val),
            ha="center",
            va=va,
            fontsize=6.8,
            color="#333333",
        )


def set_baseline_ylim(ax, vals, errs=None, baseline=0.5, pad=0.018, min_span=0.10) -> None:
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    if errs is None:
        errs = np.zeros(len(vals))
    else:
        errs = np.asarray(errs, dtype=float)
        errs = errs[~np.isnan(errs)]
    if len(vals) == 0:
        ax.set_ylim(baseline - min_span / 2, baseline + min_span / 2)
        return
    if len(errs) != len(vals):
        errs = np.zeros(len(vals))
    low = min(float(np.min(vals - errs)), baseline) - pad
    high = max(float(np.max(vals + errs)), baseline) + pad
    if high - low < min_span:
        mid = (high + low) / 2
        low = mid - min_span / 2
        high = mid + min_span / 2
    ax.set_ylim(low, high)


def plot_two_metric_context_bars(
    ax,
    data: pd.DataFrame,
    context_order: list[str],
    value_col: str,
    err_col: str | None,
    ylabel: str,
    title: str,
    value_fmt: str,
    show_error: bool = False,
) -> None:
    x = np.arange(len(context_order))
    width = 0.34
    for metric, offset in [("pearson", -width / 2), ("spearman", width / 2)]:
        sub = data[data["metric"].eq(metric)].set_index("context_label").reindex(context_order)
        vals = sub[value_col].to_numpy(dtype=float)
        errs = sub[err_col].fillna(0).to_numpy(dtype=float) if show_error and err_col else np.zeros(len(vals))
        bar_kwargs = {
            "color": COLORS[metric],
            "edgecolor": "#555555",
            "linewidth": 0.55,
            "label": METRIC_LABELS[metric],
        }
        if show_error and err_col:
            bar_kwargs.update(
                {
                    "yerr": errs,
                    "capsize": 2.4,
                    "error_kw": {"elinewidth": 0.75, "capthick": 0.75, "ecolor": "#333333"},
                }
            )
        bars = ax.bar(x + offset, vals - 0.5, width, bottom=0.5, **bar_kwargs)
        add_value_labels(ax, bars, vals, errs, fmt=value_fmt, dy=0.002)
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(context_order, rotation=28, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold", pad=13)
    ax.legend(ncol=2, loc="upper left", bbox_to_anchor=(0.0, 1.09), borderaxespad=0.0, handlelength=1.6)


def plot_five_metric_level_bars(
    ax,
    data: pd.DataFrame,
    value_col: str,
    err_col: str,
    ylabel: str,
    title: str,
    value_fmt: str,
) -> None:
    x = np.arange(len(METRIC_ORDER))
    width = 0.34
    for level, offset, color in [
        ("Tissue contexts", -width / 2, "#6B8FB3"),
        ("Adult cell types", width / 2, "#8AAE92"),
    ]:
        sub = data[data["context_level"].eq(level)].set_index("metric").reindex(METRIC_ORDER)
        vals = sub[value_col].to_numpy(dtype=float)
        errs = sub[err_col].fillna(0).to_numpy(dtype=float)
        bars = ax.bar(
            x + offset,
            vals - 0.5,
            width,
            bottom=0.5,
            yerr=errs,
            capsize=2.4,
            color=color,
            edgecolor="#555555",
            linewidth=0.55,
            error_kw={"elinewidth": 0.75, "capthick": 0.75, "ecolor": "#333333"},
            label=REVISED_CONTEXT_LABELS[level].replace("\n", " "),
        )
        add_value_labels(ax, bars, vals, errs, fmt=value_fmt, dy=0.002)
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([METRIC_LABELS[m] for m in METRIC_ORDER], rotation=22, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)


def figure_forward_revised(context_summary: pd.DataFrame, metric_summary: pd.DataFrame) -> None:
    pear_spear = context_summary[context_summary["metric"].isin(["pearson", "spearman"])].copy()
    fig = plt.figure(figsize=(9.4, 8.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.92], hspace=0.72)
    ax = fig.add_subplot(gs[0, 0])
    tissue = pear_spear[pear_spear["context_level"].eq("Tissue contexts")]
    plot_two_metric_context_bars(
        ax,
        tissue,
        TISSUE_CONTEXT_ORDER,
        "mean_auprc",
        None,
        "Mean AUPRC",
        "A. PBMC and adult tissue contexts",
        "{:.3f}",
    )
    set_baseline_ylim(ax, tissue["mean_auprc"].to_numpy(dtype=float), pad=0.018, min_span=0.10)

    ax = fig.add_subplot(gs[1, 0])
    cell = pear_spear[pear_spear["context_level"].eq("Adult cell types")]
    plot_two_metric_context_bars(
        ax,
        cell,
        CELLTYPE_CONTEXT_ORDER,
        "mean_auprc",
        None,
        "Mean AUPRC",
        "B. Adult cell-type contexts",
        "{:.3f}",
    )
    set_baseline_ylim(ax, cell["mean_auprc"].to_numpy(dtype=float), pad=0.018, min_span=0.10)

    ax = fig.add_subplot(gs[2, 0])
    plot_five_metric_level_bars(
        ax,
        metric_summary,
        "mean_auprc",
        "sd_auprc",
        "Mean AUPRC",
        "C. Five-metric comparison",
        "{:.3f}",
    )
    set_baseline_ylim(
        ax,
        metric_summary["mean_auprc"].to_numpy(dtype=float),
        metric_summary["sd_auprc"].fillna(0).to_numpy(dtype=float),
        pad=0.018,
        min_span=0.10,
    )
    fig.suptitle(
        "Curated TF–target edges retain statistical visibility across PBMC, adult tissue, and cell-type contexts",
        fontsize=11.2,
        fontweight="bold",
        y=0.995,
    )
    save_figure(fig, "integrated_figure2_forward_visibility")


def figure_reverse_revised(context_summary: pd.DataFrame, metric_summary: pd.DataFrame) -> None:
    top100_context = context_summary[context_summary["top_level"].eq("top_100")].copy()
    top100_metric = metric_summary[metric_summary["top_level"].eq("top_100")].copy()
    pear_spear = top100_context[top100_context["metric"].isin(["pearson", "spearman"])].copy()
    fig = plt.figure(figsize=(9.4, 8.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.92], hspace=0.72)

    ax = fig.add_subplot(gs[0, 0])
    tissue = pear_spear[pear_spear["context_level"].eq("Tissue contexts")]
    plot_two_metric_context_bars(
        ax,
        tissue,
        TISSUE_CONTEXT_ORDER,
        "mean_curated_positive_fraction",
        None,
        "Curated-positive fraction",
        "A. PBMC and adult tissue contexts",
        "{:.3f}",
    )
    set_baseline_ylim(
        ax,
        tissue["mean_curated_positive_fraction"].to_numpy(dtype=float),
        pad=0.020,
        min_span=0.16,
    )

    ax = fig.add_subplot(gs[1, 0])
    cell = pear_spear[pear_spear["context_level"].eq("Adult cell types")]
    plot_two_metric_context_bars(
        ax,
        cell,
        CELLTYPE_CONTEXT_ORDER,
        "mean_curated_positive_fraction",
        None,
        "Curated-positive fraction",
        "B. Adult cell-type contexts",
        "{:.3f}",
    )
    set_baseline_ylim(
        ax,
        cell["mean_curated_positive_fraction"].to_numpy(dtype=float),
        pad=0.020,
        min_span=0.16,
    )

    ax = fig.add_subplot(gs[2, 0])
    plot_five_metric_level_bars(
        ax,
        top100_metric,
        "mean_curated_positive_fraction",
        "sd_curated_positive_fraction",
        "Top100 curated-positive fraction",
        "C. Five-metric comparison",
        "{:.3f}",
    )
    set_baseline_ylim(
        ax,
        top100_metric["mean_curated_positive_fraction"].to_numpy(dtype=float),
        top100_metric["sd_curated_positive_fraction"].fillna(0).to_numpy(dtype=float),
        pad=0.020,
        min_span=0.16,
    )
    fig.suptitle(
        "Reverse consistency of high-scoring TF–gene pairs across PBMC, adult tissue, and adult cell-type contexts",
        fontsize=11.2,
        fontweight="bold",
        y=0.995,
    )
    save_figure(fig, "integrated_figure3_reverse_consistency")


def figure_hallmark_revised(hallmark: pd.DataFrame) -> None:
    top100 = hallmark[hallmark["top_level"].eq("top_100")].copy()
    fig, axes = plt.subplots(3, 1, figsize=(8.6, 8.25), sharex=False)
    for ax, group in zip(axes, HALLMARK_GROUP_ORDER):
        sub = top100[top100["context_group"].eq(group)].copy()
        sub = sub.sort_values(["support_fraction", "median_fold_enrichment"], ascending=False).head(7)
        sub = sub.sort_values("median_fold_enrichment", ascending=False)
        sub = sub.iloc[::-1]
        y = np.arange(len(sub))
        vals = sub["median_fold_enrichment"].to_numpy(dtype=float)
        bars = ax.barh(
            y,
            vals,
            color=HALLMARK_GROUP_COLORS[group],
            edgecolor="#555555",
            linewidth=0.55,
        )
        for bar, row in zip(bars, sub.itertuples(index=False)):
            label = f"{int(row.combinations_fdr_in_at_least_5_repeats)}/{int(row.combinations)}"
            ax.text(
                row.median_fold_enrichment + 0.06,
                bar.get_y() + bar.get_height() / 2,
                label,
                va="center",
                ha="left",
                fontsize=6.7,
                color="#333333",
            )
        ax.axvline(1.0, color="#333333", linestyle="--", linewidth=0.75)
        ax.set_yticks(y)
        ax.set_yticklabels([HALLMARK_DISPLAY_LABELS.get(p, p) for p in sub["pathway"]])
        ax.set_title(HALLMARK_GROUP_LABELS[group], fontweight="bold", loc="left", pad=6)
        ax.set_xlabel("Median fold enrichment")
        ax.set_xlim(0, max(3.2, float(np.nanmax(vals)) + 0.65))
    axes[0].set_ylabel("Hallmark pathway")
    fig.suptitle(
        "Hallmark enrichment among targets from high-scoring curated TF–target pairs",
        fontsize=11.0,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.01,
        0.005,
        "Bar labels show supported/analyzable combinations; support denotes FDR < 0.05 in at least 5 of 10 matched repeats.",
        fontsize=7.0,
        ha="left",
        va="bottom",
    )
    fig.subplots_adjust(top=0.925, bottom=0.08, left=0.31, right=0.98, hspace=0.55)
    save_figure(fig, "integrated_figure4_hallmark_coherence")


def figure_hallmark_context_specific_supplement() -> None:
    detail = pd.read_csv(RESULT_OUT / "revised_hallmark_analysis_units.csv", low_memory=False)
    top100 = detail[detail["top_level"].eq("top_100")].copy()
    if "context_group" not in top100.columns:
        top100["context_group"] = pd.NA
        top100.loc[top100["dataset"].isin(PBMC_DATASET_LABELS), "context_group"] = "PBMC datasets"
        top100.loc[top100["dataset"].eq("ZJU_HCL_adult_tissue"), "context_group"] = "Adult tissues"
        top100.loc[top100["dataset"].eq("ZJU_HCL_adult_celltype"), "context_group"] = "Adult cell types"
    top100 = top100[top100["context_group"].notna()].copy()
    top100["supported"] = top100["repeats_fdr_lt_0_05"] >= 5
    context = (
        top100.groupby(["context_group", "context_label", "pathway"], dropna=False)
        .agg(
            combinations=("supported", "size"),
            supported=("supported", "sum"),
            support_fraction=("supported", "mean"),
            median_fold_enrichment=("median_fold_enrichment", "median"),
        )
        .reset_index()
    )
    plot_df = context.copy()
    plot_df["context_display"] = plot_df["context_label"].str.replace("\n", " ", regex=False)
    panels = [
        ("PBMC datasets", label.replace("\n", " "))
        for label in TISSUE_CONTEXT_ORDER[:4]
    ] + [
        ("Adult tissues", label.replace("\n", " "))
        for label in TISSUE_CONTEXT_ORDER[4:]
    ] + [
        ("Adult cell types", label)
        for label in CELLTYPE_CONTEXT_ORDER
    ]

    fig = plt.figure(figsize=(11.5, 13.0))
    gs = fig.add_gridspec(7, 3, hspace=0.82, wspace=0.85)
    row_slots = [
        (0, 0), (0, 1), (0, 2), (1, 0),
        (2, 0), (2, 1), (2, 2), (3, 0),
        (4, 0), (4, 1), (4, 2), (5, 0), (5, 1), (5, 2), (6, 0),
    ]
    axes = []
    for (group, context_label), (row, col) in zip(panels, row_slots):
        ax = fig.add_subplot(gs[row, col])
        axes.append(ax)
        sub = plot_df[
            plot_df["context_group"].eq(group)
            & plot_df["context_display"].eq(context_label)
        ].copy()
        sub = sub.sort_values(["support_fraction", "median_fold_enrichment"], ascending=False).head(4)
        sub = sub.sort_values("median_fold_enrichment", ascending=True)
        y = np.arange(len(sub))
        vals = sub["median_fold_enrichment"].to_numpy(dtype=float)
        ax.barh(
            y,
            vals,
            color=HALLMARK_GROUP_COLORS[group],
            edgecolor="#555555",
            linewidth=0.5,
        )
        for yi, row_data in enumerate(sub.itertuples(index=False)):
            label = f"{int(row_data.supported)}/{int(row_data.combinations)}"
            ax.text(
                float(row_data.median_fold_enrichment) + 0.05,
                yi,
                label,
                va="center",
                ha="left",
                fontsize=5.6,
                color="#333333",
            )
        ax.axvline(1.0, color="#333333", linestyle="--", linewidth=0.65)
        ax.set_yticks(y)
        ax.set_yticklabels([HALLMARK_DISPLAY_LABELS.get(p, p) for p in sub["pathway"]], fontsize=5.8)
        ax.tick_params(axis="x", labelsize=5.8)
        ax.set_title(context_label, fontweight="bold", fontsize=7.2, pad=3.5)
        if row == 6 or (row == 5 and col > 0):
            ax.set_xlabel("Median FE", fontsize=6.2)
        xmax = max(2.2, float(np.nanmax(vals)) + 0.75) if len(vals) else 2.2
        ax.set_xlim(0, xmax)

    used_slots = set(row_slots[: len(panels)])
    for row in range(7):
        for col in range(3):
            if (row, col) not in used_slots:
                ax = fig.add_subplot(gs[row, col])
                ax.axis("off")

    fig.text(0.02, 0.947, "PBMC datasets", fontsize=8.2, fontweight="bold", color="#333333")
    fig.text(0.02, 0.675, "Adult tissues", fontsize=8.2, fontweight="bold", color="#333333")
    fig.text(0.02, 0.403, "Adult cell types", fontsize=8.2, fontweight="bold", color="#333333")
    fig.suptitle(
        "Context-specific Hallmark enrichment among targets from high-scoring curated TF–target pairs",
        fontsize=10.8,
        fontweight="bold",
        y=0.992,
    )
    fig.text(
        0.02,
        0.012,
        "Each panel shows the top four pathways for that context. Bar labels show supported/analyzable combinations; FE = fold enrichment.",
        fontsize=7.0,
        ha="left",
        va="bottom",
    )
    fig.subplots_adjust(top=0.925, bottom=0.045, left=0.19, right=0.985)
    save_figure(fig, "integrated_figure4_context_specific_supplement")


def main() -> None:
    setup()
    forward = load_forward_tables()
    forward_context_summary, forward_metric_summary = summarize_forward_revised(forward)
    reverse = load_reverse_tables()
    reverse_context_summary, reverse_metric_summary = summarize_reverse_revised(reverse)
    hallmark = summarize_hallmark_revised()
    write_supplementary_workbook(
        forward_metric_summary,
        forward_context_summary,
        reverse_metric_summary,
        reverse_context_summary,
        hallmark,
    )
    figure_forward_revised(forward_context_summary, forward_metric_summary)
    figure_reverse_revised(reverse_context_summary, reverse_metric_summary)
    figure_hallmark_revised(hallmark)
    print(f"[write] {RESULT_OUT}")
    print(f"[write] {FIG_OUT}")


if __name__ == "__main__":
    main()
