# -*- coding: utf-8 -*-
"""
Create final result bar plots for the repositioned GRN visibility manuscript.

Figures:
- Figure 3: Forward validation main result across datasets and metrics.
- Figure 4: Reverse validation known-positive fraction among high-scoring pairs.
- Figure 5: Pathway enrichment of high-scoring positive targets.

The script reads only current refined result tables and writes new final figures.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
FORWARD = ROOT / "results" / "forward_validation_refined"
REVERSE = ROOT / "results" / "reverse_validation_from_refined_pairs"
OUT = ROOT / "figures" / "final"

DATASET_ORDER = [
    "PBMC10k",
    "Kang_IFN_beta",
    "GSE178429_IFN_gamma",
    "NYGC_multimodal_PBMC",
    "GSE126030_T_cells",
]
DATASET_LABELS = {
    "PBMC10k": "PBMC10k",
    "Kang_IFN_beta": "Kang\nIFN-beta",
    "GSE178429_IFN_gamma": "GSE178429\nIFN-gamma",
    "NYGC_multimodal_PBMC": "NYGC\nPBMC",
    "GSE126030_T_cells": "GSE126030\nT cells",
}
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
COLORS = {
    "pearson": "#5b7fa3",
    "spearman": "#7aa88f",
    "mutual_information": "#a9b4bf",
    "coexpression_probability": "#c9b27f",
    "codetection_odds_ratio": "#9f7f9f",
    "top_100": "#5b7fa3",
    "top_500": "#9bb6a5",
    "ci": "#8c8c8c",
    "Hallmark": "#83AA8D",
    "Reactome": "#7D9DB8",
}


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        path = OUT / f"{name}.{ext}"
        if ext == "png":
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
        print(path)
    plt.close(fig)


def ci95(values: pd.Series) -> float:
    values = values.dropna()
    if len(values) <= 1:
        return 0.0
    return 1.96 * values.std(ddof=1) / np.sqrt(len(values))


def add_bar_labels(ax, bars, values, errors=None, dy=0.002, fmt="{:.3f}") -> None:
    if errors is None:
        errors = np.zeros(len(values), dtype=float)
    for bar, value, err in zip(bars, values, errors):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + float(err) + dy,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=7.4,
            color="#333333",
            rotation=0,
        )


def figure3_forward_main(forward_repeat: pd.DataFrame, forward_summary: pd.DataFrame) -> None:
    df = forward_repeat[forward_repeat["metric"].isin(["pearson", "spearman"])].copy()
    summary = (
        df.groupby(["dataset", "metric"], dropna=False)
        .agg(mean_auprc=("auprc", "mean"), ci=("auprc", ci95), n=("auprc", "size"))
        .reset_index()
    )

    x = np.arange(len(DATASET_ORDER))
    width = 0.34
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax = axes[0]
    offsets = {"pearson": -width / 2, "spearman": width / 2}
    for metric in ["pearson", "spearman"]:
        sub = summary[summary["metric"].eq(metric)].set_index("dataset").reindex(DATASET_ORDER).reset_index()
        vals = sub["mean_auprc"].to_numpy()
        errs = sub["ci"].fillna(0).to_numpy()
        bars = ax.bar(
            x + offsets[metric],
            vals,
            width,
            yerr=errs,
            capsize=3.0,
            error_kw={"elinewidth": 0.85, "capthick": 0.85, "ecolor": "#333333"},
            color=COLORS[metric],
            edgecolor="#555555",
            linewidth=0.7,
            label=METRIC_LABELS[metric],
        )
        add_bar_labels(ax, bars, vals, errors=errs, dy=0.002)

    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([DATASET_LABELS[d] for d in DATASET_ORDER])
    ax.set_ylabel("AUPRC")
    ax.set_title("A. Dataset-level visibility")
    ax.set_ylim(0.48, max(0.60, summary["mean_auprc"].max() + 0.025))
    ax.legend(title="", frameon=False, ncol=2, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    metric_df = forward_summary.copy()
    metric_summary = (
        metric_df.groupby("metric", dropna=False)
        .agg(mean_auprc=("auprc_mean", "mean"), ci=("auprc_mean", ci95), n_units=("auprc_mean", "size"))
        .reindex(METRIC_ORDER)
        .reset_index()
    )
    ax = axes[1]
    mx = np.arange(len(METRIC_ORDER))
    vals = metric_summary["mean_auprc"].to_numpy()
    errs = metric_summary["ci"].fillna(0).to_numpy()
    bars = ax.bar(
        mx,
        vals,
        yerr=errs,
        capsize=3.0,
        error_kw={"elinewidth": 0.85, "capthick": 0.85, "ecolor": "#333333"},
        color=[COLORS[m] for m in METRIC_ORDER],
        edgecolor="#555555",
        linewidth=0.7,
    )
    add_bar_labels(ax, bars, vals, errors=errs, dy=0.002)
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=0.9)
    ax.set_xticks(mx)
    ax.set_xticklabels([METRIC_LABELS[m] for m in METRIC_ORDER], rotation=22, ha="right")
    ax.set_ylabel("Mean AUPRC")
    ax.set_title("B. Five-metric comparison")
    ax.set_ylim(0.48, max(0.58, metric_summary["mean_auprc"].max() + 0.025))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle(
        "Curated TF-target edges show weak but reproducible statistical visibility",
        y=1.02,
        fontsize=11.5,
        fontweight="bold",
    )
    save(fig, "figure3_forward_main_barplot")


def figure4_reverse_known_fraction(reverse_summary: pd.DataFrame) -> None:
    df = reverse_summary[reverse_summary["top_level"].isin(["top_100", "top_500"])].copy()
    summary = (
        df.groupby(["metric", "top_level"], dropna=False)
        .agg(mean_fraction=("known_positive_fraction_mean", "mean"), ci=("known_positive_fraction_mean", ci95), n_units=("known_positive_fraction_mean", "size"))
        .reset_index()
    )
    x = np.arange(len(METRIC_ORDER))
    width = 0.34
    offsets = {"top_100": -width / 2, "top_500": width / 2}
    labels = {"top_100": "Top 100", "top_500": "Top 500"}

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for top_level in ["top_100", "top_500"]:
        sub = summary[summary["top_level"].eq(top_level)].set_index("metric").reindex(METRIC_ORDER).reset_index()
        vals = sub["mean_fraction"].to_numpy()
        errs = sub["ci"].fillna(0).to_numpy()
        bars = ax.bar(
            x + offsets[top_level],
            vals,
            width,
            yerr=errs,
            capsize=3.0,
            error_kw={"elinewidth": 0.85, "capthick": 0.85, "ecolor": "#333333"},
            color=COLORS[top_level],
            edgecolor="#555555",
            linewidth=0.7,
            label=labels[top_level],
        )
        add_bar_labels(ax, bars, vals, errors=errs, dy=0.002, fmt="{:.2f}")

    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([METRIC_LABELS[m] for m in METRIC_ORDER], rotation=18, ha="right")
    ax.set_ylabel("Known-positive fraction")
    ax.set_title("Top-ranked pairs are enriched for known regulatory edges")
    ax.set_ylim(0.46, max(0.68, summary["mean_fraction"].max() + 0.04))
    ax.legend(title="", frameon=False, ncol=2, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "figure4_reverse_known_fraction_barplot")


def figure5_pathway_enrichment(pathway_summary: pd.DataFrame) -> None:
    plot = pathway_summary[
        pathway_summary["top_level"].eq("top_100")
        & pathway_summary["edge_set"].eq("trrust_dorothea_intersection")
    ].copy()
    pathway_order = [
        ("Hallmark", "IL-6/JAK/STAT3 Signaling", "IL-6/JAK/STAT3 signaling"),
        ("Hallmark", "Interferon Gamma Response", "Interferon gamma response"),
        ("Hallmark", "Interferon Alpha Response", "Interferon alpha response"),
        ("Hallmark", "Inflammatory Response", "Inflammatory response"),
        ("Reactome", "Interferon Alpha Beta Signaling", "Interferon alpha/beta signaling"),
        ("Reactome", "Interferon Signaling", "Interferon signaling"),
        ("Hallmark", "IL-2/STAT5 Signaling", "IL-2/STAT5 signaling"),
    ]
    rows = []
    for collection, pathway, label in pathway_order:
        sub = plot[plot["collection"].eq(collection) & plot["pathway"].eq(pathway)]
        if sub.empty:
            continue
        rows.append(
            {
                "collection": collection,
                "label": label,
                "median_fe": sub["median_fold_enrichment"].median(),
                "consistent_units": int((sub["repeats_fdr_lt_0_05"] >= 5).sum()),
                "total_units": len(sub),
            }
        )
    path_df = pd.DataFrame(rows).sort_values("median_fe", ascending=True)

    fig, ax = plt.subplots(figsize=(7.4, 4.85))
    y = np.arange(len(path_df))
    ax.barh(
        y,
        path_df["median_fe"],
        color=[COLORS[c] for c in path_df["collection"]],
        edgecolor="#424242",
        linewidth=0.55,
        height=0.58,
    )
    ax.axvline(1.0, color="#666666", linewidth=0.9, linestyle="--")
    for yi, row in zip(y, path_df.itertuples(index=False)):
        ax.text(row.median_fe + 0.07, yi, f"{row.consistent_units}/{row.total_units}", va="center", ha="left", fontsize=8.5, color="#222222")
    ax.set_yticks(y)
    ax.set_yticklabels(path_df["label"], fontsize=8.8)
    ax.set_xlabel("Median fold enrichment of target genes", fontsize=9.4)
    ax.set_xlim(0, max(3.05, path_df["median_fe"].max() + 0.75))
    ax.tick_params(axis="x", labelsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#444444")
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_title("Pathway enrichment of high-scoring positive targets", fontsize=10.8, fontweight="bold", pad=14)
    ax.text(0, 1.012, "Top 100 ranked database-supported positive pairs; five datasets", transform=ax.transAxes, fontsize=7.9, color="#555555", ha="left", va="bottom")
    legend = [
        Patch(facecolor=COLORS["Hallmark"], edgecolor="#424242", label="Hallmark"),
        Patch(facecolor=COLORS["Reactome"], edgecolor="#424242", label="Reactome"),
    ]
    ax.legend(handles=legend, frameon=False, loc="lower right", fontsize=8.2, bbox_to_anchor=(1.0, 0.05))
    fig.subplots_adjust(left=0.335, right=0.985, top=0.79, bottom=0.16)
    save(fig, "figure5_reverse_pathway_context")


def main() -> None:
    setup()
    forward_repeat = pd.read_csv(FORWARD / "forward_refined_by_repeat.csv")
    forward_summary = pd.read_csv(FORWARD / "forward_refined_full_summary.csv")
    reverse_summary = pd.read_csv(REVERSE / "reverse_refined_known_overlap_summary.csv")
    pathway_summary = pd.read_csv(REVERSE / "reverse_refined_pathway_validation_summary.csv")
    figure3_forward_main(forward_repeat, forward_summary)
    figure4_reverse_known_fraction(reverse_summary)
    figure5_pathway_enrichment(pathway_summary)


if __name__ == "__main__":
    main()
