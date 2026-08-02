"""
Formal undirected association visibility analysis.

Pearson and Spearman are scored as absolute correlations. Mutual information,
co-expression probability, and co-detection odds ratio keep their existing
definitions. Background matching keeps the TF fixed, matches target expression
and detection bins, and excludes all known targets from TRRUST union DoRothEA
A/B.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "formal_abs_association"
TABLES = ROOT / "results" / "tables"
N_REPEATS = 10
N_BINS = 8
SEED = 20260602
DETECTION_THRESHOLD = 0.05
T_95_DF9 = 2.2621571627409915

FORMAL_METRICS = [
    "pearson",
    "spearman",
    "mutual_information",
    "coexpression_probability",
    "codetection_odds_ratio",
]


PBMC_EDGE_SETS = {
    "trrust": ROOT / "data" / "ground_truth" / "trrust_edges_detection_ge_05pct.csv",
    "dorothea_ab": ROOT / "data" / "ground_truth" / "dorothea_ab_edges_detection_ge_05pct.csv",
    "trrust_dorothea_intersection": ROOT / "data" / "ground_truth" / "trrust_dorothea_ab_intersection_detection_ge_05pct.csv",
    "trrust_dorothea_union": ROOT / "data" / "ground_truth" / "trrust_dorothea_ab_union_detection_ge_05pct.csv",
}

INTERSECTION_BASE = ROOT / "data" / "ground_truth" / "trrust_dorothea_ab_intersection_edges_standardized.csv"
KNOWN_FILES = [
    ROOT / "data" / "ground_truth" / "trrust_edges_standardized.csv",
    ROOT / "data" / "ground_truth" / "dorothea_ab_edges_standardized.csv",
]

DATASETS = [
    {
        "dataset": "PBMC10k",
        "condition": "all",
        "h5ad": ROOT / "data" / "processed" / "pbmc10k_normalized.h5ad",
        "gene_summary": TABLES / "pbmc10k_gene_summary.csv",
        "edge_sets": PBMC_EDGE_SETS,
    },
    {
        "dataset": "Kang_IFN_beta",
        "condition": "ctrl",
        "h5ad": ROOT / "data" / "processed" / "kang2018_ctrl_normalized.h5ad",
        "gene_summary": TABLES / "kang2018_ctrl_gene_summary.csv",
        "edge_sets": {"trrust_dorothea_intersection": INTERSECTION_BASE},
    },
    {
        "dataset": "Kang_IFN_beta",
        "condition": "stim",
        "h5ad": ROOT / "data" / "processed" / "kang2018_stim_normalized.h5ad",
        "gene_summary": TABLES / "kang2018_stim_gene_summary.csv",
        "edge_sets": {"trrust_dorothea_intersection": INTERSECTION_BASE},
    },
    {
        "dataset": "GSE178429_IFN_gamma",
        "condition": "ctrl_6h",
        "h5ad": ROOT / "data" / "processed" / "gse178429_ctrl_6h_normalized.h5ad",
        "gene_summary": TABLES / "gse178429_ctrl_6h_gene_summary.csv",
        "edge_sets": {"trrust_dorothea_intersection": INTERSECTION_BASE},
    },
    {
        "dataset": "GSE178429_IFN_gamma",
        "condition": "ifng_6h",
        "h5ad": ROOT / "data" / "processed" / "gse178429_ifng_6h_normalized.h5ad",
        "gene_summary": TABLES / "gse178429_ifng_6h_gene_summary.csv",
        "edge_sets": {"trrust_dorothea_intersection": INTERSECTION_BASE},
    },
]


def standardize_mode(mode: object) -> str:
    value = str(mode).strip().lower()
    if value == "activation":
        return "Activation"
    if value == "repression":
        return "Repression"
    return "Unknown"


def standardize_edges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out["target"] = out["target"].astype(str).str.upper()
    if "mode" not in out.columns:
        out["mode"] = "Unknown"
    out["mode"] = out["mode"].map(standardize_mode)
    return out


def add_bins(df: pd.DataFrame, col: str) -> pd.Series:
    ranks = df[col].rank(method="first")
    return pd.qcut(ranks, q=min(N_BINS, len(df)), labels=False, duplicates="drop")


def standardized_mean_difference(x: pd.Series, y: pd.Series) -> float:
    vx = x.var(ddof=1)
    vy = y.var(ddof=1)
    pooled = np.sqrt((vx + vy) / 2)
    if pooled == 0 or np.isnan(pooled):
        return 0.0
    return float((x.mean() - y.mean()) / pooled)


def load_known_targets() -> dict[str, set[str]]:
    edges = pd.concat([standardize_edges(pd.read_csv(p))[["tf", "target"]] for p in KNOWN_FILES], ignore_index=True)
    edges = edges.drop_duplicates(["tf", "target"])
    return edges.groupby("tf")["target"].apply(set).to_dict()


def load_gene_summary(path: Path) -> pd.DataFrame:
    genes = pd.read_csv(path)
    genes["gene"] = genes["gene"].astype(str).str.upper()
    return genes.drop_duplicates("gene").reset_index(drop=True)


def prepare_positives(path: Path, genes: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
    raw = standardize_edges(pd.read_csv(path))
    before = len(raw)
    stats = genes.set_index("gene")
    present = set(stats.index)
    pos = raw[raw["tf"].isin(present) & raw["target"].isin(present)].copy()

    for side, col in [("tf", "tf"), ("target", "target")]:
        pos[f"{side}_detection_rate"] = pos[col].map(stats["detection_rate"])
        pos[f"{side}_mean_counts"] = pos[col].map(stats["mean_counts"])

    pos = pos[
        (pos["tf_detection_rate"] >= DETECTION_THRESHOLD)
        & (pos["target_detection_rate"] >= DETECTION_THRESHOLD)
    ].drop_duplicates(["tf", "target"], keep="first")
    return pos.reset_index(drop=True), before, len(pos)


def generate_negatives(
    positives: pd.DataFrame,
    genes: pd.DataFrame,
    known_targets: dict[str, set[str]],
    dataset: str,
    condition: str,
    edge_set: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    genes = genes.copy()
    genes["expr_bin"] = add_bins(genes, "mean_counts")
    genes["detect_bin"] = add_bins(genes, "detection_rate")
    stats = genes.set_index("gene")
    rng = np.random.default_rng(SEED)
    rows = []
    failures = 0

    for repeat in range(N_REPEATS):
        for edge_id, row in positives.reset_index(drop=True).iterrows():
            tf = row["tf"]
            target = row["target"]
            if target not in stats.index:
                failures += 1
                continue
            ts = stats.loc[target]
            candidates = genes[
                (genes["expr_bin"] == ts["expr_bin"])
                & (genes["detect_bin"] == ts["detect_bin"])
            ]["gene"].tolist()
            excluded = known_targets.get(tf, set()) | {tf, target}
            candidates = [g for g in candidates if g not in excluded]

            relaxed = 0
            if not candidates:
                relaxed = 1
                candidates = genes[
                    (genes["expr_bin"].between(ts["expr_bin"] - 1, ts["expr_bin"] + 1))
                    & (genes["detect_bin"].between(ts["detect_bin"] - 1, ts["detect_bin"] + 1))
                ]["gene"].tolist()
                candidates = [g for g in candidates if g not in excluded]

            if not candidates:
                failures += 1
                continue

            neg_target = rng.choice(candidates)
            ns = stats.loc[neg_target]
            rows.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "repeat": repeat,
                    "edge_id": edge_id,
                    "tf": tf,
                    "positive_target": target,
                    "target": neg_target,
                    "mode": row.get("mode", "Unknown"),
                    "relaxed_bin": relaxed,
                    "pos_target_mean_counts": row["target_mean_counts"],
                    "neg_target_mean_counts": ns["mean_counts"],
                    "pos_target_detection_rate": row["target_detection_rate"],
                    "neg_target_detection_rate": ns["detection_rate"],
                    "tf_mean_counts": row["tf_mean_counts"],
                    "tf_detection_rate": row["tf_detection_rate"],
                }
            )

    neg = pd.DataFrame(rows)
    pos_pairs = set(map(tuple, positives[["tf", "target"]].to_numpy()))
    neg_pairs = set(map(tuple, neg[["tf", "target"]].to_numpy())) if len(neg) else set()
    known_pairs = {(tf, target) for tf, targets in known_targets.items() for target in targets}
    qc = {
        "dataset": dataset,
        "condition": condition,
        "edge_set": edge_set,
        "positive_edges_after_filter_before_dedup": int(len(positives)),
        "positive_edges_after_dedup": int(len(positives)),
        "n_negative_edges": int(len(neg)),
        "repeats": N_REPEATS,
        "per_repeat_positive_edges": int(len(positives)),
        "per_repeat_negative_edges_min": int(neg.groupby("repeat").size().min()) if len(neg) else 0,
        "per_repeat_negative_edges_max": int(neg.groupby("repeat").size().max()) if len(neg) else 0,
        "failed_matches": int(failures),
        "relaxed_fraction": float(neg["relaxed_bin"].mean()) if len(neg) else np.nan,
        "smd_target_mean_counts": standardized_mean_difference(neg["pos_target_mean_counts"], neg["neg_target_mean_counts"]) if len(neg) else np.nan,
        "smd_target_detection_rate": standardized_mean_difference(neg["pos_target_detection_rate"], neg["neg_target_detection_rate"]) if len(neg) else np.nan,
        "positive_negative_overlap_unique": len(pos_pairs & neg_pairs),
        "negative_overlap_with_known_union_unique": len(neg_pairs & known_pairs),
        "negative_unique_edges": len(neg_pairs),
        "negative_duplicate_rows": int(len(neg) - len(neg_pairs)),
    }
    return neg, qc


def matching_qc_from_tables(
    positives: pd.DataFrame,
    negatives: pd.DataFrame,
    known_targets: dict[str, set[str]],
    dataset: str,
    condition: str,
    edge_set: str,
) -> dict[str, object]:
    pos_pairs = set(map(tuple, positives[["tf", "target"]].to_numpy()))
    neg_pairs = set(map(tuple, negatives[["tf", "target"]].to_numpy())) if len(negatives) else set()
    known_pairs = {(tf, target) for tf, targets in known_targets.items() for target in targets}
    return {
        "dataset": dataset,
        "condition": condition,
        "edge_set": edge_set,
        "positive_edges_after_filter_before_dedup": int(len(positives)),
        "positive_edges_after_dedup": int(len(positives)),
        "n_negative_edges": int(len(negatives)),
        "repeats": int(negatives["repeat"].nunique()) if len(negatives) else N_REPEATS,
        "per_repeat_positive_edges": int(len(positives)),
        "per_repeat_negative_edges_min": int(negatives.groupby("repeat").size().min()) if len(negatives) else 0,
        "per_repeat_negative_edges_max": int(negatives.groupby("repeat").size().max()) if len(negatives) else 0,
        "failed_matches": int(len(positives) * N_REPEATS - len(negatives)),
        "relaxed_fraction": float(negatives["relaxed_bin"].mean()) if len(negatives) else np.nan,
        "smd_target_mean_counts": standardized_mean_difference(negatives["pos_target_mean_counts"], negatives["neg_target_mean_counts"]) if len(negatives) else np.nan,
        "smd_target_detection_rate": standardized_mean_difference(negatives["pos_target_detection_rate"], negatives["neg_target_detection_rate"]) if len(negatives) else np.nan,
        "positive_negative_overlap_unique": len(pos_pairs & neg_pairs),
        "negative_overlap_with_known_union_unique": len(neg_pairs & known_pairs),
        "negative_unique_edges": len(neg_pairs),
        "negative_duplicate_rows": int(len(negatives) - len(neg_pairs)),
    }


def decode_names(values: np.ndarray) -> list[str]:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def read_csr(group: h5py.Group) -> sparse.csr_matrix:
    shape = tuple(group.attrs["shape"])
    return sparse.csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)


def read_var_names(handle: h5py.File) -> list[str]:
    index_name = handle["var"].attrs.get("_index", "gene")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return decode_names(handle["var"][index_name][:])


def read_h5ad_submatrix(h5ad: Path, needed_genes: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    with h5py.File(h5ad, "r") as handle:
        var_names = read_var_names(handle)
        gene_index = {str(g).upper(): i for i, g in enumerate(var_names)}
        available = [g for g in needed_genes if g in gene_index]
        original_indices = [gene_index[g] for g in available]
        x = read_csr(handle["X"])[:, original_indices].toarray().astype(np.float32, copy=False)
        counts = read_csr(handle["layers"]["counts"])[:, original_indices].toarray() > 0
    return x, counts, available


def mi_from_bins(x: np.ndarray, y: np.ndarray) -> float:
    counts = np.bincount(
        x.astype(np.int16) * N_BINS + y.astype(np.int16),
        minlength=N_BINS * N_BINS,
    ).astype(np.float64)
    pxy = counts.reshape(N_BINS, N_BINS)
    pxy /= pxy.sum()
    px = pxy.sum(axis=1, keepdims=True)
    py = pxy.sum(axis=0, keepdims=True)
    expected = px @ py
    mask = pxy > 0
    return float(np.sum(pxy[mask] * np.log(pxy[mask] / expected[mask])))


def compute_pair_metrics(h5ad: Path, pairs: pd.DataFrame) -> pd.DataFrame:
    unique = standardize_edges(pairs[["tf", "target"]]).drop_duplicates(["tf", "target"]).reset_index(drop=True)
    needed = sorted(set(unique["tf"]) | set(unique["target"]))
    x_norm, x_binary, available = read_h5ad_submatrix(h5ad, needed)
    sub_index = {gene: idx for idx, gene in enumerate(available)}
    valid = unique[unique["tf"].isin(sub_index) & unique["target"].isin(sub_index)].reset_index(drop=True)
    left = valid["tf"].map(sub_index).to_numpy()
    right = valid["target"].map(sub_index).to_numpy()
    n_cells = x_norm.shape[0]

    mean = x_norm.mean(axis=0)
    std = x_norm.std(axis=0)
    detected = x_binary.sum(axis=0)

    ranks = np.empty(x_norm.shape, dtype=np.float32)
    bins = np.empty(x_norm.shape, dtype=np.uint8)
    for start in range(0, x_norm.shape[1], 500):
        end = min(start + 500, x_norm.shape[1])
        r = rankdata(x_norm[:, start:end], axis=0)
        ranks[:, start:end] = r.astype(np.float32, copy=False)
        bins[:, start:end] = np.floor((r - 1) * N_BINS / n_cells).clip(0, N_BINS - 1).astype(np.uint8)

    rank_mean = ranks.mean(axis=0)
    rank_std = ranks.std(axis=0)

    out = valid.copy()
    for metric in FORMAL_METRICS:
        out[metric] = np.nan

    for start in range(0, len(valid), 2000):
        end = min(start + 2000, len(valid))
        li = left[start:end]
        ri = right[start:end]

        mean_xy = (x_norm[:, li] * x_norm[:, ri]).mean(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out.loc[start:end - 1, "pearson"] = (mean_xy - mean[li] * mean[ri]) / (std[li] * std[ri])

        rank_mean_xy = (ranks[:, li] * ranks[:, ri]).mean(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out.loc[start:end - 1, "spearman"] = (rank_mean_xy - rank_mean[li] * rank_mean[ri]) / (rank_std[li] * rank_std[ri])

        both = np.logical_and(x_binary[:, li], x_binary[:, ri]).sum(axis=0).astype(np.float32)
        x_only = detected[li].astype(np.float32) - both
        y_only = detected[ri].astype(np.float32) - both
        neither = n_cells - both - x_only - y_only
        out.loc[start:end - 1, "coexpression_probability"] = both / n_cells
        out.loc[start:end - 1, "codetection_odds_ratio"] = (
            (both + 0.5) * (neither + 0.5)
        ) / ((x_only + 0.5) * (y_only + 0.5))

    mi = np.full(len(valid), np.nan, dtype=np.float32)
    for idx, (li, ri) in enumerate(zip(left, right)):
        mi[idx] = mi_from_bins(bins[:, li], bins[:, ri])
    out["mutual_information"] = mi
    return out


def score_values(df: pd.DataFrame, metric: str, score_definition: str) -> pd.Series:
    if score_definition == "formal_abs" and metric in {"pearson", "spearman"}:
        return df[metric].abs()
    return df[metric]


def cliffs_delta(pos: pd.Series, neg: pd.Series) -> float:
    x = pos.to_numpy()
    y = neg.to_numpy()
    ranks = rankdata(np.concatenate([x, y]))
    n_pos = len(x)
    n_neg = len(y)
    u = ranks[:n_pos].sum() - n_pos * (n_pos + 1) / 2
    return float((2 * u / (n_pos * n_neg)) - 1)


def enrichment_at_5pct(values: pd.Series, labels: pd.Series) -> float:
    k = max(1, int(np.ceil(len(values) * 0.05)))
    top = values.sort_values(ascending=False).index[:k]
    return float(labels.loc[top].mean() / labels.mean())


def evaluate(
    dataset: str,
    condition: str,
    edge_set: str,
    positives: pd.DataFrame,
    negatives: pd.DataFrame,
    pair_metrics: pd.DataFrame,
    score_definition: str,
) -> pd.DataFrame:
    rows = []
    for repeat in sorted(negatives["repeat"].unique()):
        pos = positives.reset_index(names="edge_id")[["edge_id", "tf", "target", "mode"]].copy()
        pos["repeat"] = repeat
        pos["label"] = 1
        neg = negatives[negatives["repeat"].eq(repeat)][["repeat", "edge_id", "tf", "target", "mode"]].copy()
        neg["label"] = 0
        test = pd.concat([pos, neg], ignore_index=True)
        test = test.merge(pair_metrics, on=["tf", "target"], how="left")
        for metric in FORMAL_METRICS:
            sub = test[["label", metric]].replace([np.inf, -np.inf], np.nan).dropna()
            if sub.empty or sub["label"].nunique() < 2:
                continue
            scores = score_values(sub, metric, score_definition)
            y = sub["label"].astype(int)
            pos_scores = scores.loc[y == 1]
            neg_scores = scores.loc[y == 0]
            rows.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "score_definition": score_definition,
                    "metric": metric,
                    "repeat": int(repeat),
                    "n_positive": int((y == 1).sum()),
                    "n_negative": int((y == 0).sum()),
                    "random_auprc_baseline": float(y.mean()),
                    "auroc": roc_auc_score(y, scores),
                    "auprc": average_precision_score(y, scores),
                    "enrichment_at_5pct": enrichment_at_5pct(scores, y),
                    "positive_median": float(pos_scores.median()),
                    "negative_median": float(neg_scores.median()),
                    "median_difference": float(pos_scores.median() - neg_scores.median()),
                    "cliffs_delta": cliffs_delta(pos_scores, neg_scores),
                }
            )
    return pd.DataFrame(rows)


def add_ci(mean: pd.Series, sd: pd.Series, n: pd.Series) -> tuple[pd.Series, pd.Series]:
    half_width = T_95_DF9 * sd / np.sqrt(n)
    return mean - half_width, mean + half_width


def summarize(by_repeat: pd.DataFrame) -> pd.DataFrame:
    summary = (
        by_repeat.groupby(["dataset", "condition", "edge_set", "score_definition", "metric"])
        .agg(
            n_repeats=("repeat", "nunique"),
            n_positive_mean=("n_positive", "mean"),
            n_negative_mean=("n_negative", "mean"),
            random_auprc_baseline=("random_auprc_baseline", "mean"),
            auroc_mean=("auroc", "mean"),
            auroc_sd=("auroc", "std"),
            auprc_mean=("auprc", "mean"),
            auprc_sd=("auprc", "std"),
            enrichment_at_5pct_mean=("enrichment_at_5pct", "mean"),
            enrichment_at_5pct_sd=("enrichment_at_5pct", "std"),
            median_difference_mean=("median_difference", "mean"),
            median_difference_sd=("median_difference", "std"),
            cliffs_delta_mean=("cliffs_delta", "mean"),
            cliffs_delta_sd=("cliffs_delta", "std"),
            positive_median_mean=("positive_median", "mean"),
            negative_median_mean=("negative_median", "mean"),
        )
        .reset_index()
    )
    for col in ["auroc", "auprc", "enrichment_at_5pct", "median_difference", "cliffs_delta"]:
        lo, hi = add_ci(summary[f"{col}_mean"], summary[f"{col}_sd"], summary["n_repeats"])
        summary[f"{col}_ci95_low"] = lo
        summary[f"{col}_ci95_high"] = hi
    return summary


def main() -> None:
    """Explain why this helper module is not a public analysis entry point."""
    print(
        "\n".join(
            [
                "42_formal_abs_association_main_results.py is a helper module.",
                "It provides shared functions for metric calculation and evaluation.",
                "The manuscript-facing refined analyses are generated or summarized by scripts 43-46.",
                "Running the original full matching workflow requires raw or intermediate expression matrices that are not included in this GitHub package.",
            ]
        )
    )


if __name__ == "__main__":
    main()
