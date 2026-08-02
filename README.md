# TF-target statistical visibility in single-cell expression data

This repository contains the analysis code and processed result tables for a two-part evaluation of curated TF-target regulatory relationships in single-cell RNA-seq expression matrices.

## Study scope

1. Forward validation: curated TF-target edges from TRRUST v2 and DoRothEA A/B are compared with strictly matched background pairs.
2. Reverse consistency: high-scoring TF-gene pairs from the evaluated pair space are checked for curated-positive fraction and pathway coherence.

This project does not perform de novo GRN reconstruction, causal validation, TF-activity modelling, machine-learning classification, or large-language-model-based analysis.

## Official result sources

- `results/forward_validation_refined/`
- `results/reverse_validation_from_refined_pairs/`
- `results/formal_abs_association_refined_matching/`
- `results/server_added_refined_matching/`

Pearson and Spearman refer to absolute correlations in the formal analyses.

The formal summary includes five single-cell data sources: PBMC10k, NYGC multimodal PBMC, Kang IFN-beta PBMC, GSE178429 IFN-gamma PBMC, and GSE126030 T-cell tissue contexts.

## Main scripts

- `scripts/42_formal_abs_association_main_results.py`: helper functions for pair-level metric calculation and evaluation. It is not a public full-analysis entry point.
- `scripts/43_known_edge_visibility_refined_matching.py`: refined matched-background workflow for the original local datasets. Full regeneration requires raw or intermediate expression matrices that are not included in this GitHub package.
- `scripts/44_forward_validation_refined_summary.py`: rebuilds the combined forward-validation summary from included result tables.
- `scripts/45_reverse_validation_from_refined_pairs.py`: rebuilds reverse-consistency summaries from included refined pair results. This script can also generate a large repeat-level top-pair detail file that is intentionally not included in this package.
- `scripts/46_reverse_refined_pathway_validation.py`: checks or rebuilds pathway-validation summaries from included reverse-validation outputs.
- `scripts/51_make_final_result_barplots.py`: creates optional local barplots from included summary tables.

## Quick checks

After installing the dependencies in `requirements.txt`, the following commands should run without requiring raw expression matrices:

```bash
python scripts/42_formal_abs_association_main_results.py
python scripts/44_forward_validation_refined_summary.py
python scripts/46_reverse_refined_pathway_validation.py
python scripts/51_make_final_result_barplots.py
```

The reverse-consistency script can be rerun when detailed top-pair outputs are needed:

```bash
python scripts/45_reverse_validation_from_refined_pairs.py
```

This command may take longer than the quick checks and may create `results/reverse_validation_from_refined_pairs/reverse_refined_top_pairs_by_repeat.csv`, which is excluded from version control.

## Data and file-size notes

Raw and processed single-cell expression matrices are not included because they are large and publicly available from the original sources cited in the manuscript. The repository includes processed result tables needed to verify the manuscript numbers.

The detailed file `results/reverse_validation_from_refined_pairs/reverse_refined_top_pairs_by_repeat.csv` is excluded from this GitHub package because it is a large repeat-level detail file. It is not required to verify the manuscript summary tables, because the included reverse-consistency summary files contain the values reported in the paper and supplement.

## Manuscript boundary

The manuscript claim is limited to statistical visibility and reverse consistency of curated TF-target edges in single-cell expression space. The analysis does not claim causal regulation, direct TF binding, confirmed non-regulatory background edges, or de novo regulatory-network reconstruction.
