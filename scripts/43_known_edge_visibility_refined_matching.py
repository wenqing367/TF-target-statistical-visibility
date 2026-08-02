"""
Official refined known-edge visibility pipeline.

This script regenerates the current strict matched-background analysis used for
the project question: whether known TF-target edges show stronger undirected
statistical association than matched background TF-gene pairs in single-cell
expression matrices.

It intentionally keeps the project scope narrow:
- Pearson and Spearman are scored as absolute correlations.
- Mutual information, co-expression probability, and co-detection odds ratio
  reuse the definitions from 42_formal_abs_association_main_results.py.
- Background edges keep the same TF as the positive edge.
- Background targets are matched on target mean expression and detection rate.
- Known targets from TRRUST union DoRothEA A/B, the TF itself, and the original
  target are excluded.
- Each repeat uses unique negative TF-target pairs.
- Only positives matched in all 10 repeats are retained as the strict common
  support subset.

The script will not overwrite an existing output directory unless --overwrite
is supplied. This protects results/formal_abs_association_refined_matching/
from accidental changes.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HELPER_PATH = ROOT / "scripts" / "42_formal_abs_association_main_results.py"
DEFAULT_OUT = ROOT / "results" / "formal_abs_association_refined_matching"

CALIPER_Z = 0.25
TOP_K_NEAREST = 5
SEED = 20260602


spec = importlib.util.spec_from_file_location("formal_abs_helpers", HELPER_PATH)
helpers = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(helpers)


def add_standardized_background(genes: pd.DataFrame) -> pd.DataFrame:
    out = genes.copy()
    for col, zcol in [("mean_counts", "z_mean_counts"), ("detection_rate", "z_detection_rate")]:
        sd = out[col].std(ddof=0)
        out[zcol] = 0.0 if sd == 0 or np.isnan(sd) else (out[col] - out[col].mean()) / sd
    return out


def candidate_lists(
    positives: pd.DataFrame,
    genes: pd.DataFrame,
    known_targets: dict[str, set[str]],
) -> tuple[dict[int, pd.DataFrame], dict[int, str]]:
    gene_names = genes["gene"].to_numpy()
    z_mean = genes["z_mean_counts"].to_numpy()
    z_detect = genes["z_detection_rate"].to_numpy()
    lists: dict[int, pd.DataFrame] = {}
    reasons: dict[int, str] = {}

    for edge_id, row in positives.reset_index(drop=True).iterrows():
        tf = row["tf"]
        target = row["target"]
        excluded = known_targets.get(tf, set()) | {tf, target}
        dz_mean = z_mean - row["target_z_mean_counts"]
        dz_detect = z_detect - row["target_z_detection_rate"]
        mask = (np.abs(dz_mean) <= CALIPER_Z) & (np.abs(dz_detect) <= CALIPER_Z)
        if not np.any(mask):
            reasons[edge_id] = "no_gene_within_caliper"
            lists[edge_id] = pd.DataFrame(columns=["gene", "distance"])
            continue

        cand = pd.DataFrame(
            {
                "gene": gene_names[mask],
                "distance": np.sqrt(dz_mean[mask] ** 2 + dz_detect[mask] ** 2),
            }
        )
        cand = cand[~cand["gene"].isin(excluded)].sort_values(["distance", "gene"]).reset_index(drop=True)
        if cand.empty:
            reasons[edge_id] = "all_caliper_candidates_excluded"
        lists[edge_id] = cand
    return lists, reasons


def refined_match(
    positives: pd.DataFrame,
    genes: pd.DataFrame,
    known_targets: dict[str, set[str]],
    dataset: str,
    condition: str,
    edge_set: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    positives = positives.copy().reset_index(drop=True)
    genes = add_standardized_background(genes)
    stats = genes.set_index("gene")
    positives["target_z_mean_counts"] = positives["target"].map(stats["z_mean_counts"])
    positives["target_z_detection_rate"] = positives["target"].map(stats["z_detection_rate"])

    lists, base_reasons = candidate_lists(positives, genes, known_targets)
    rng = np.random.default_rng(SEED)
    order_base = sorted(
        range(len(positives)),
        key=lambda i: (len(lists[i]), positives.loc[i, "tf"], positives.loc[i, "target"]),
    )

    rows = []
    failure_records = []
    for repeat in range(helpers.N_REPEATS):
        used_pairs: set[tuple[str, str]] = set()
        ordered = []
        for count in sorted({len(lists[i]) for i in order_base}):
            block = [i for i in order_base if len(lists[i]) == count]
            rng.shuffle(block)
            ordered.extend(block)

        for edge_id in ordered:
            row = positives.loc[edge_id]
            tf = row["tf"]
            cand = lists[edge_id]
            if cand.empty:
                failure_records.append(
                    {
                        "dataset": dataset,
                        "condition": condition,
                        "edge_set": edge_set,
                        "repeat": repeat,
                        "edge_id": edge_id,
                        "tf": tf,
                        "target": row["target"],
                        "reason": base_reasons.get(edge_id, "no_candidate"),
                    }
                )
                continue

            available = cand[~cand["gene"].map(lambda g: (tf, g) in used_pairs)]
            if available.empty:
                failure_records.append(
                    {
                        "dataset": dataset,
                        "condition": condition,
                        "edge_set": edge_set,
                        "repeat": repeat,
                        "edge_id": edge_id,
                        "tf": tf,
                        "target": row["target"],
                        "reason": "caliper_candidates_exhausted_by_unique_constraint",
                    }
                )
                continue

            nearest = available.head(TOP_K_NEAREST)
            weights = np.exp(-nearest["distance"].to_numpy() / 0.05)
            weights = weights / weights.sum()
            picked = nearest.loc[rng.choice(nearest.index.to_numpy(), p=weights)]
            neg_target = str(picked["gene"])
            used_pairs.add((tf, neg_target))
            ns = stats.loc[neg_target]
            rows.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "edge_set": edge_set,
                    "repeat": repeat,
                    "edge_id": edge_id,
                    "tf": tf,
                    "positive_target": row["target"],
                    "target": neg_target,
                    "mode": row.get("mode", "Unknown"),
                    "match_distance": float(picked["distance"]),
                    "caliper_z": CALIPER_Z,
                    "relaxed_bin": 0,
                    "pos_target_mean_counts": row["target_mean_counts"],
                    "neg_target_mean_counts": ns["mean_counts"],
                    "pos_target_detection_rate": row["target_detection_rate"],
                    "neg_target_detection_rate": ns["detection_rate"],
                    "tf_mean_counts": row["tf_mean_counts"],
                    "tf_detection_rate": row["tf_detection_rate"],
                }
            )

    negatives_initial = pd.DataFrame(rows)
    matched_counts = (
        negatives_initial.groupby("edge_id")["repeat"].nunique()
        if len(negatives_initial)
        else pd.Series(dtype=int)
    )
    common_edge_ids = set(matched_counts[matched_counts == helpers.N_REPEATS].index.astype(int))
    positives_common = positives.loc[sorted(common_edge_ids)].copy().reset_index(drop=True)
    edge_id_map = {old: new for new, old in enumerate(sorted(common_edge_ids))}

    negatives = negatives_initial[negatives_initial["edge_id"].isin(common_edge_ids)].copy()
    negatives["edge_id_original"] = negatives["edge_id"]
    negatives["edge_id"] = negatives["edge_id"].map(edge_id_map).astype(int)
    positives_common["edge_id_original"] = sorted(common_edge_ids)

    failures = pd.DataFrame(failure_records)
    if len(failures):
        failures["in_common_supported_subset"] = failures["edge_id"].isin(common_edge_ids)
    return positives_common, negatives.reset_index(drop=True), failures


def qc_tables(
    positives_original_n: int,
    positives: pd.DataFrame,
    negatives: pd.DataFrame,
    failures: pd.DataFrame,
    known_targets: dict[str, set[str]],
    dataset: str,
    condition: str,
    edge_set: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    known_pairs = {(tf, target) for tf, targets in known_targets.items() for target in targets}
    pos_pairs = set(map(tuple, positives[["tf", "target"]].to_numpy()))
    neg_pairs = set(map(tuple, negatives[["tf", "target"]].to_numpy())) if len(negatives) else set()
    per_rows = []

    for repeat, sub in negatives.groupby("repeat"):
        internal_dup = len(sub) - len(sub[["tf", "target"]].drop_duplicates())
        per_rows.append(
            {
                "dataset": dataset,
                "condition": condition,
                "edge_set": edge_set,
                "repeat": int(repeat),
                "positive_edges": len(positives),
                "negative_edges": len(sub),
                "internal_duplicate_tf_target": internal_dup,
                "failed_matches": len(positives) - len(sub),
                "relaxed_fraction": float(sub["relaxed_bin"].mean()) if len(sub) else np.nan,
                "positive_negative_overlap_unique": len(pos_pairs & set(map(tuple, sub[["tf", "target"]].to_numpy()))),
                "negative_overlap_with_known_union_unique": len(set(map(tuple, sub[["tf", "target"]].to_numpy())) & known_pairs),
                "smd_target_mean_counts": helpers.standardized_mean_difference(
                    sub["pos_target_mean_counts"], sub["neg_target_mean_counts"]
                )
                if len(sub)
                else np.nan,
                "smd_target_detection_rate": helpers.standardized_mean_difference(
                    sub["pos_target_detection_rate"], sub["neg_target_detection_rate"]
                )
                if len(sub)
                else np.nan,
            }
        )

    per_repeat = pd.DataFrame(per_rows)
    reason_counts = (
        failures.drop_duplicates(["edge_id", "reason"])["reason"].value_counts().to_dict()
        if len(failures)
        else {}
    )
    summary = {
        "dataset": dataset,
        "condition": condition,
        "edge_set": edge_set,
        "caliper_z": CALIPER_Z,
        "original_positive_edges": positives_original_n,
        "matched_positive_edges_common": len(positives),
        "matching_coverage": len(positives) / positives_original_n if positives_original_n else np.nan,
        "unmatched_positive_edges": positives_original_n - len(positives),
        "unmatched_reasons": "; ".join(f"{k}:{v}" for k, v in reason_counts.items()),
        "repeats": int(per_repeat["repeat"].nunique()) if len(per_repeat) else 0,
        "per_repeat_positive_edges": len(positives),
        "per_repeat_negative_min": int(per_repeat["negative_edges"].min()) if len(per_repeat) else 0,
        "per_repeat_negative_max": int(per_repeat["negative_edges"].max()) if len(per_repeat) else 0,
        "internal_duplicate_total": int(per_repeat["internal_duplicate_tf_target"].sum()) if len(per_repeat) else 0,
        "internal_duplicate_max": int(per_repeat["internal_duplicate_tf_target"].max()) if len(per_repeat) else 0,
        "failed_matches_total_final": int(per_repeat["failed_matches"].sum()) if len(per_repeat) else 0,
        "relaxed_fraction_max": float(per_repeat["relaxed_fraction"].max()) if len(per_repeat) else np.nan,
        "positive_negative_overlap_unique": len(pos_pairs & neg_pairs),
        "negative_overlap_with_known_union_unique": len(neg_pairs & known_pairs),
        "smd_target_mean_counts_mean": float(per_repeat["smd_target_mean_counts"].mean()) if len(per_repeat) else np.nan,
        "smd_target_mean_counts_max_abs": float(per_repeat["smd_target_mean_counts"].abs().max()) if len(per_repeat) else np.nan,
        "smd_target_detection_rate_mean": float(per_repeat["smd_target_detection_rate"].mean()) if len(per_repeat) else np.nan,
        "smd_target_detection_rate_max_abs": float(per_repeat["smd_target_detection_rate"].abs().max()) if len(per_repeat) else np.nan,
    }
    return per_repeat, summary


def run_analysis(out: Path, use_existing_components: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    known_targets = helpers.load_known_targets()
    all_qc_rows = []
    all_per_repeat_qc = []
    all_failures = []
    eval_rows = []

    for config in helpers.DATASETS:
        dataset = config["dataset"]
        condition = config["condition"]
        genes = helpers.load_gene_summary(config["gene_summary"])
        analyses = []
        all_pairs = []

        for edge_set, edge_path in config["edge_sets"].items():
            positives_raw, _before, after = helpers.prepare_positives(edge_path, genes)
            prefix = f"{dataset}_{condition}_{edge_set}".replace("/", "_")
            pos_path = out / f"{prefix}_refined_positives.csv"
            neg_path = out / f"{prefix}_refined_negatives.csv"
            fail_path = out / f"{prefix}_refined_unmatched_failures.csv"

            if use_existing_components and pos_path.exists() and neg_path.exists() and fail_path.exists():
                positives = pd.read_csv(pos_path)
                negatives = pd.read_csv(neg_path)
                failures = pd.read_csv(fail_path)
            else:
                positives, negatives, failures = refined_match(
                    positives_raw, genes, known_targets, dataset, condition, edge_set
                )
                positives.to_csv(pos_path, index=False)
                negatives.to_csv(neg_path, index=False)
                failures.to_csv(fail_path, index=False)

            per_qc, summary_qc = qc_tables(
                after, positives, negatives, failures, known_targets, dataset, condition, edge_set
            )
            all_qc_rows.append(summary_qc)
            all_per_repeat_qc.append(per_qc)
            if len(failures):
                all_failures.append(failures)
            analyses.append((edge_set, positives, negatives))
            all_pairs.append(pd.concat([positives[["tf", "target"]], negatives[["tf", "target"]]], ignore_index=True))
            print(f"[matched] {prefix}: original={after} common={len(positives)} negatives={len(negatives)}")

        pair_input = pd.concat(all_pairs, ignore_index=True).drop_duplicates(["tf", "target"])
        pair_prefix = f"{dataset}_{condition}".replace("/", "_")
        pair_path = out / f"{pair_prefix}_refined_pair_metrics.csv"
        if use_existing_components and pair_path.exists():
            pair_metrics = pd.read_csv(pair_path)
        else:
            pair_metrics = helpers.compute_pair_metrics(config["h5ad"], pair_input)
            pair_metrics.to_csv(pair_path, index=False)

        for edge_set, positives, negatives in analyses:
            eval_rows.append(
                helpers.evaluate(dataset, condition, edge_set, positives, negatives, pair_metrics, "formal_abs")
            )

    by_repeat = pd.concat(eval_rows, ignore_index=True)
    summary = helpers.summarize(by_repeat)
    summary = summary[summary["score_definition"].eq("formal_abs")].reset_index(drop=True)
    qc = pd.DataFrame(all_qc_rows)
    per_qc = pd.concat(all_per_repeat_qc, ignore_index=True)
    failures_all = pd.concat(all_failures, ignore_index=True) if all_failures else pd.DataFrame()

    by_repeat.to_csv(out / "refined_formal_abs_by_repeat.csv", index=False)
    summary.to_csv(out / "refined_formal_abs_summary.csv", index=False)
    qc.to_csv(out / "refined_matching_qc_summary.csv", index=False)
    per_qc.to_csv(out / "refined_matching_qc_by_repeat.csv", index=False)
    failures_all.to_csv(out / "refined_unmatched_failures_all.csv", index=False)
    print(f"[write] {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUT,
        help="Output directory. Defaults to results/formal_abs_association_refined_matching/.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow writing into a non-empty output directory.",
    )
    parser.add_argument(
        "--reuse-existing-components",
        action="store_true",
        help="Reuse existing refined positive/negative/pair-metric CSV files in the output directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out = args.output_dir
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise SystemExit(
            f"{out} already contains files. Use --overwrite to refresh it, "
            "or pass --output-dir to write an independent reproducibility run."
        )
    run_analysis(out, use_existing_components=args.reuse_existing_components)


if __name__ == "__main__":
    main()
