# Reverse Validation From Refined Pairs

Scope: rank only the TF-gene pairs already evaluated in the strict refined forward validation.
No new candidate universe is generated and no expression metric is recomputed.

- Analysis units: 12
- Top-pair rows: 524800
- Overlap summary rows: 180
- Pathway summary rows: 1162

Balanced random baseline: each repeat contains 1:1 curated positives and matched background pairs, so the expected curated-positive fraction among top-ranked pairs is 0.5 if the metric carries no ranking signal.

Pearson/Spearman top curated-positive fractions:
             dataset condition                     edge_set   metric top_level  known_positive_fraction_mean  delta_vs_balanced_baseline_mean  repeats_above_0_5
 GSE178429_IFN_gamma   ctrl_6h trrust_dorothea_intersection  pearson   top_100                        0.6420                           0.1420                 10
 GSE178429_IFN_gamma   ctrl_6h trrust_dorothea_intersection  pearson   top_500                        0.5152                           0.0152                 10
 GSE178429_IFN_gamma   ctrl_6h trrust_dorothea_intersection spearman   top_100                        0.6530                           0.1530                 10
 GSE178429_IFN_gamma   ctrl_6h trrust_dorothea_intersection spearman   top_500                        0.5122                           0.0122                 10
 GSE178429_IFN_gamma   ifng_6h trrust_dorothea_intersection  pearson   top_100                        0.5880                           0.0880                 10
 GSE178429_IFN_gamma   ifng_6h trrust_dorothea_intersection  pearson   top_500                        0.5502                           0.0502                 10
 GSE178429_IFN_gamma   ifng_6h trrust_dorothea_intersection spearman   top_100                        0.6230                           0.1230                 10
 GSE178429_IFN_gamma   ifng_6h trrust_dorothea_intersection spearman   top_500                        0.5652                           0.0652                 10
       Kang_IFN_beta      ctrl trrust_dorothea_intersection  pearson   top_100                        0.5380                           0.0380                 10
       Kang_IFN_beta      ctrl trrust_dorothea_intersection  pearson   top_500                        0.5052                           0.0052                 10
       Kang_IFN_beta      ctrl trrust_dorothea_intersection spearman   top_100                        0.5590                           0.0590                 10
       Kang_IFN_beta      ctrl trrust_dorothea_intersection spearman   top_500                        0.5084                           0.0084                  9
       Kang_IFN_beta      stim trrust_dorothea_intersection  pearson   top_100                        0.5470                           0.0470                 10
       Kang_IFN_beta      stim trrust_dorothea_intersection  pearson   top_500                        0.4980                          -0.0020                  3
       Kang_IFN_beta      stim trrust_dorothea_intersection spearman   top_100                        0.5230                           0.0230                  8
       Kang_IFN_beta      stim trrust_dorothea_intersection spearman   top_500                        0.5064                           0.0064                  9
NYGC_multimodal_PBMC       all                  dorothea_ab  pearson   top_100                        0.5780                           0.0780                 10
NYGC_multimodal_PBMC       all                  dorothea_ab  pearson   top_500                        0.5176                           0.0176                 10
NYGC_multimodal_PBMC       all                  dorothea_ab spearman   top_100                        0.5610                           0.0610                 10
NYGC_multimodal_PBMC       all                  dorothea_ab spearman   top_500                        0.5164                           0.0164                 10
NYGC_multimodal_PBMC       all                       trrust  pearson   top_100                        0.6350                           0.1350                 10
NYGC_multimodal_PBMC       all                       trrust  pearson   top_500                        0.5452                           0.0452                 10
NYGC_multimodal_PBMC       all                       trrust spearman   top_100                        0.6410                           0.1410                 10
NYGC_multimodal_PBMC       all                       trrust spearman   top_500                        0.5210                           0.0210                 10
NYGC_multimodal_PBMC       all trrust_dorothea_intersection  pearson   top_100                        0.5500                           0.0500                 10
NYGC_multimodal_PBMC       all trrust_dorothea_intersection  pearson   top_500                        0.5254                           0.0254                 10
NYGC_multimodal_PBMC       all trrust_dorothea_intersection spearman   top_100                        0.5510                           0.0510                 10
NYGC_multimodal_PBMC       all trrust_dorothea_intersection spearman   top_500                        0.5126                           0.0126                  9
NYGC_multimodal_PBMC       all        trrust_dorothea_union  pearson   top_100                        0.6170                           0.1170                 10
NYGC_multimodal_PBMC       all        trrust_dorothea_union  pearson   top_500                        0.5322                           0.0322                 10
NYGC_multimodal_PBMC       all        trrust_dorothea_union spearman   top_100                        0.5830                           0.0830                 10
NYGC_multimodal_PBMC       all        trrust_dorothea_union spearman   top_500                        0.5224                           0.0224                 10
             PBMC10k       all                  dorothea_ab  pearson   top_100                        0.5900                           0.0900                 10
             PBMC10k       all                  dorothea_ab  pearson   top_500                        0.5150                           0.0150                 10
             PBMC10k       all                  dorothea_ab spearman   top_100                        0.5900                           0.0900                 10
             PBMC10k       all                  dorothea_ab spearman   top_500                        0.5052                           0.0052                  9
             PBMC10k       all                       trrust  pearson   top_100                        0.6330                           0.1330                 10
             PBMC10k       all                       trrust  pearson   top_500                        0.5234                           0.0234                 10
             PBMC10k       all                       trrust spearman   top_100                        0.6490                           0.1490                 10
             PBMC10k       all                       trrust spearman   top_500                        0.5160                           0.0160                  9
             PBMC10k       all trrust_dorothea_intersection  pearson   top_100                        0.5540                           0.0540                 10
             PBMC10k       all trrust_dorothea_intersection  pearson   top_500                        0.5198                           0.0198                 10
             PBMC10k       all trrust_dorothea_intersection spearman   top_100                        0.5370                           0.0370                 10
             PBMC10k       all trrust_dorothea_intersection spearman   top_500                        0.5180                           0.0180                 10
             PBMC10k       all        trrust_dorothea_union  pearson   top_100                        0.6270                           0.1270                 10
             PBMC10k       all        trrust_dorothea_union  pearson   top_500                        0.5194                           0.0194                 10
             PBMC10k       all        trrust_dorothea_union spearman   top_100                        0.6280                           0.1280                 10
             PBMC10k       all        trrust_dorothea_union spearman   top_500                        0.5194                           0.0194                 10
