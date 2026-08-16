import pandas as pd
from _paths import PROJECT_ROOT as ROOT

IN=ROOT/'results'/'zju_cellatlas_adult_tissue_refined_matching_by_context'
OUT=ROOT/'results'/'extension_forward_visibility'/'hcl_adult_tissue'
OUT.mkdir(parents=True, exist_ok=True)
files={
 'refined_formal_abs_summary.csv':'zju_hcl_adult_tissue_refined_formal_abs_summary_combined.csv',
 'refined_formal_abs_by_repeat.csv':'zju_hcl_adult_tissue_refined_formal_abs_by_repeat_combined.csv',
 'refined_matching_qc_summary.csv':'zju_hcl_adult_tissue_refined_matching_qc_summary_combined.csv',
 'refined_matching_qc_by_repeat.csv':'zju_hcl_adult_tissue_refined_matching_qc_by_repeat_combined.csv',
}
for fname,outname in files.items():
    parts=[]
    for p in sorted(IN.glob(f'*/{fname}')):
        df=pd.read_csv(p)
        df.insert(0,'run_dir',p.parent.name)
        parts.append(df)
    if not parts:
        raise SystemExit(f'No files found for {fname}')
    all_df=pd.concat(parts, ignore_index=True)
    all_df.to_csv(OUT/outname, index=False)
    print('[write]', OUT/outname, all_df.shape)
summary=pd.read_csv(OUT/'zju_hcl_adult_tissue_refined_formal_abs_summary_combined.csv')
metric_col='cliffs_delta_mean' if 'cliffs_delta_mean' in summary.columns else 'cliffs_delta'
over=summary.groupby('metric').agg(
    units=('auprc_mean','size'),
    auprc_mean=('auprc_mean','mean'),
    auprc_median=('auprc_mean','median'),
    auprc_gt_05=('auprc_mean', lambda x:int((x>0.5).sum())),
    auroc_mean=('auroc_mean','mean'),
    auroc_median=('auroc_mean','median'),
    auroc_gt_05=('auroc_mean', lambda x:int((x>0.5).sum())),
    cliffs_delta_mean=(metric_col,'mean'),
).reset_index()
over.to_csv(OUT/'zju_hcl_adult_tissue_metric_overview.csv', index=False)
print('[write]', OUT/'zju_hcl_adult_tissue_metric_overview.csv', over.shape)
qc=pd.read_csv(OUT/'zju_hcl_adult_tissue_refined_matching_qc_summary_combined.csv')
qc_over=pd.DataFrame([{ 
    'qc_rows':len(qc),
    'matching_coverage_min':qc['matching_coverage'].min(),
    'matching_coverage_max':qc['matching_coverage'].max(),
    'internal_duplicate_total_max':qc['internal_duplicate_total'].max(),
    'internal_duplicate_max_max':qc['internal_duplicate_max'].max(),
    'failed_matches_total_final_max':qc['failed_matches_total_final'].max(),
    'relaxed_fraction_max':qc['relaxed_fraction_max'].max(),
    'positive_negative_overlap_unique_max':qc['positive_negative_overlap_unique'].max(),
    'negative_overlap_with_known_union_unique_max':qc['negative_overlap_with_known_union_unique'].max(),
    'smd_target_mean_counts_max_abs':qc['smd_target_mean_counts_max_abs'].max(),
    'smd_target_detection_rate_max_abs':qc['smd_target_detection_rate_max_abs'].max(),
}])
qc_over.to_csv(OUT/'zju_hcl_adult_tissue_qc_overview.csv', index=False)
print('[write]', OUT/'zju_hcl_adult_tissue_qc_overview.csv', qc_over.shape)
print(over.to_string(index=False))
print(qc_over.to_string(index=False))
