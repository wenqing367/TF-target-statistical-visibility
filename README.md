# TF-target statistical visibility in single-cell expression data

This repository contains processed result tables and reproducible plotting code for evaluating the statistical visibility of curated TF-target regulatory relationships in single-cell expression matrices.

The included manuscript-level workflow was verified with Python 3.8.8 and the pinned package versions in `requirements.txt`.

## Study scope

The current manuscript-level analysis has two parts.

1. Forward visibility: curated TF-target edges from TRRUST v2 and DoRothEA A/B are compared with strictly matched background TF-gene pairs.
2. Reverse consistency and pathway coherence: high-scoring TF-gene pairs are checked for curated-positive fraction, and target genes from high-scoring curated TF-target pairs are tested for Hallmark pathway enrichment.

The analysis is limited to expression-level statistical visibility. It does not perform de novo GRN reconstruction, causal validation, direct TF-binding inference, TF-activity modelling, machine-learning classification, or large-language-model-based analysis.

## Current manuscript contexts

The integrated manuscript results use three context groups.

- PBMC datasets: PBMC10k, Kang IFN-beta PBMC, GSE178429 IFN-gamma PBMC, and NYGC multimodal PBMC.
- Adult tissues from HCL: adult colon, adult kidney, adult liver, and adult lung.
- Adult cell types from HCL: AT2 cell, endothelial cell, enterocyte, epithelial cell, fibroblast, Loop of Henle, and smooth muscle cell.

## Main result sources

The current manuscript-level source tables are in:

- `results/integrated_paper_outputs/`

The integrated tables are generated from:

- `results/forward_validation_refined/`
- `results/reverse_validation_from_refined_pairs/`
- `results/formal_abs_association_refined_matching/`
- `results/server_added_refined_matching/nygc_multimodal_pbmc_refined_matching/`
- `results/extension_forward_visibility/`
- `results/extension_reverse_consistency/`
- `results/extension_hallmark_coherence/`

Pearson and Spearman refer to absolute correlations in the formal analyses.

The bottom-level analysis unit is `condition × edge-set`. Edge sets and matched-background repeats are not biological replicates. Kang and GSE178429 conditions are analyzed separately before being summarized to their dataset-level final contexts.

Reverse consistency is evaluated only within the balanced matched evaluation sets; it is not an unrestricted TF-by-gene search. Fixed Top100 is the main threshold. Strict Top500 and top-5% are sensitivity analyses. A strict fixed threshold is not analyzable when fewer than the requested number of ranked pairs is available; such units are recorded as `NA`, not as full-take substitutes.

## Figures

Current integrated figures are written to:

- `figures/integrated/`

The final manuscript figure files included in this release are:

- `integrated_figure1_study_framework_with_text.*`
- `integrated_figure2_forward_visibility.*`
- `integrated_figure3_reverse_consistency.*`
- `integrated_figure4_hallmark_coherence.*`

Figures 2-4 were renumbered from the integrated plotting outputs to match the current manuscript. PNG files are the raster versions used for the manuscript to avoid SVG font and label displacement in Microsoft Word. SVG files keep editable text.

Script 52 regenerates the manuscript-level summary tables and the data-derived Figures 2-4.

## Scripts

- `scripts/42_formal_abs_association_main_results.py`: shared helper functions for pair-level metric calculation and evaluation.
- `scripts/43_known_edge_visibility_refined_matching.py`: refined matched-background workflow for the original local datasets. Full regeneration requires raw or intermediate expression matrices that are not included in this GitHub package.
- `scripts/44_forward_validation_refined_summary.py`: rebuilds the original forward-validation summary from included refined result tables.
- `scripts/45_reverse_validation_from_refined_pairs.py`: rebuilds reverse-consistency summaries from included refined pair results. It can also generate a large repeat-level top-pair detail file that is intentionally excluded from version control.
- `scripts/46_reverse_refined_pathway_validation.py`: checks or rebuilds pathway-validation summaries from included reverse-validation outputs.
- `scripts/52_make_integrated_extension_tables_figures.py`: current manuscript-level script. It combines the PBMC, adult tissue, and adult cell-type result tables and regenerates integrated summary tables and Figures 2-4.
- `scripts/53_sanitize_release_paths.py`: replaces machine-specific paths in the three release inventory CSV files with repository-relative paths.
- `scripts/extension/`: scripts used to prepare and summarize adult HCL tissue, adult HCL cell-type, extension reverse-consistency, and Hallmark-wide pathway-coherence analyses. Full matched-pair regeneration requires the corresponding raw or intermediate single-cell matrices, which are not included here.

Scripts 45-46 retain an earlier PBMC-only Hallmark/Reactome verification workflow for auditability. Reactome outputs from those scripts are not used in the current manuscript, Figure 4, or Supplementary Table S5; the manuscript-level pathway result is Hallmark-only and is generated by script 52.

## Installation

Create an isolated Python 3.8 environment and install the pinned manuscript-level dependencies:

```bash
python -m venv .venv
# Windows
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -r requirements.txt
# Linux/macOS
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Optional full regeneration from HCL `.h5ad` inputs also requires:

```bash
python -m pip install -r requirements-full.txt
```

## Quick checks

After installing dependencies in `requirements.txt`, these commands should run without raw expression matrices:

```bash
python scripts/42_formal_abs_association_main_results.py
python scripts/44_forward_validation_refined_summary.py
python scripts/46_reverse_refined_pathway_validation.py
python scripts/52_make_integrated_extension_tables_figures.py
python scripts/53_sanitize_release_paths.py
```

Expected manuscript-level outputs from script 52 include 75 forward context-metric rows, 10 forward group-metric rows, 225 reverse context-metric-threshold rows, 30 reverse group-metric-threshold rows, and the compact `integrated_supplementary_tables.xlsx` workbook with sheets S1-S5. Strict Top500 is valid in 38 of 56 analysis units and 14 of 15 final biological contexts; adult colon has no analyzable strict Top500 edge set.

The reverse-consistency script can be rerun when detailed top-pair outputs are needed:

```bash
python scripts/45_reverse_validation_from_refined_pairs.py
```

This command may take longer and may create `results/reverse_validation_from_refined_pairs/reverse_refined_top_pairs_by_repeat.csv`, which is excluded from version control.

## Optional full regeneration

The raw HCL matrices and some upstream matching inputs are not distributed in this repository. The extension scripts no longer contain machine-specific absolute paths. Configure external inputs with environment variables before using them:

```bash
export GRN_PROJECT_ROOT=/path/to/this/repository
export GRN_MATCHING_ROOT=/path/to/upstream/matching/workspace
export GRN_HALLMARK_GMT=/path/to/h.all.v2023.1.Hs.symbols.gmt
```

On Windows PowerShell, use `$env:GRN_PROJECT_ROOT`, `$env:GRN_MATCHING_ROOT`, and `$env:GRN_HALLMARK_GMT`. `GRN_MATCHING_ROOT` must contain `scripts/43_known_edge_visibility_refined_matching.py` and the four detection-filtered edge-set files expected by the two `run_zju_*_refined.py` scripts. Inventory paths written by the preparation and reverse-consistency scripts are repository-relative.

The `h5ad`, `gene_summary`, and `metric_file` entries in the inventory CSV files are provenance paths for upstream artifacts that are not distributed in this processed-results release. They are not required for the quick checks above. To rerun the upstream HCL preparation or matching stages, generate or supply those artifacts at the recorded relative paths first.

## Data and file-size notes

Raw and processed single-cell expression matrices are not included because they are large and publicly available from the original sources cited in the manuscript. This repository includes processed result tables needed to verify the manuscript numbers and regenerate integrated figures.

Large repeat-level top-pair detail files are not included. The included summary, repeat-level, quality-control, and Hallmark summary tables are sufficient for checking the reported integrated results.

No included file exceeds GitHub's 100 MB per-file limit. Before creating a public release, run the quick checks, verify that `git status` contains only intended changes, and confirm that the repository URL in the manuscript is publicly accessible.

## Interpretation boundary

Matched background pairs are comparator pairs for statistical evaluation, not verified non-regulatory edges. The manuscript claim is limited to statistical visibility and reverse consistency of curated TF-target edges in single-cell expression space.
