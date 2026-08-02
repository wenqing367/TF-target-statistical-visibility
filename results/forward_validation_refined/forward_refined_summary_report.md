# Forward Validation Refined Summary

Scope: known TF-target edges versus strict expression/detection matched background edges.
Pearson and Spearman are absolute correlations. No matching or metric recomputation is performed here.

- Analysis units: 28
- Units passing completion/QC checks: 28/28
- Dataset-condition groups: 10
- Metric result rows: 140
- Metric row counts: {'codetection_odds_ratio': 28, 'coexpression_probability': 28, 'mutual_information': 28, 'pearson': 28, 'spearman': 28}
- Absolute Pearson/Spearman rows with AUPRC > 0.5: 56/56

Important dataset note: GSE126030 is represented by tissue-specific T-cell h5ad files. The server files do not contain condition/stim/sample metadata, so these results are tissue validation results, not anti-CD3/CD28 stimulation comparisons.

Output tables:
- forward_refined_key_results.csv
- forward_refined_qc_summary.csv
- forward_refined_by_repeat.csv
- forward_refined_qc_by_repeat.csv
- forward_refined_completion_check.csv