# Reverse Validation From Refined Pairs

Scope: rank only the TF-gene pairs already evaluated in the strict refined forward validation.
No new candidate universe is generated and no expression metric is recomputed.

- Analysis units: 28
- Top-pair rows: 1141300
- Overlap summary rows: 420
- Pathway summary rows: 2758

Balanced random baseline: each repeat contains 1:1 curated positives and matched background pairs, so the expected curated-positive fraction among top-ranked pairs is 0.5 if the metric carries no ranking signal.

Pearson/Spearman top curated-positive fractions:
             dataset   condition                     edge_set   metric top_level  known_positive_fraction_mean  delta_vs_balanced_baseline_mean  repeats_above_0_5
   GSE126030_T_cells       blood                  dorothea_ab  pearson   top_100                        0.6770                           0.1770                 10
   GSE126030_T_cells       blood                  dorothea_ab  pearson   top_500                        0.6120                           0.1120                 10
   GSE126030_T_cells       blood                  dorothea_ab spearman   top_100                        0.5800                           0.0800                 10
   GSE126030_T_cells       blood                  dorothea_ab spearman   top_500                        0.5700                           0.0700                 10
   GSE126030_T_cells       blood                       trrust  pearson   top_100                        0.5910                           0.0910                 10
   GSE126030_T_cells       blood                       trrust  pearson   top_500                        0.5412                           0.0412                 10
   GSE126030_T_cells       blood                       trrust spearman   top_100                        0.5590                           0.0590                 10
   GSE126030_T_cells       blood                       trrust spearman   top_500                        0.5018                           0.0018                  5
   GSE126030_T_cells       blood trrust_dorothea_intersection  pearson   top_100                        0.5910                           0.0910                 10
   GSE126030_T_cells       blood trrust_dorothea_intersection  pearson   top_500                        0.5122                           0.0122                 10
   GSE126030_T_cells       blood trrust_dorothea_intersection spearman   top_100                        0.5230                           0.0230                  9
   GSE126030_T_cells       blood trrust_dorothea_intersection spearman   top_500                        0.5072                           0.0072                  9
   GSE126030_T_cells       blood        trrust_dorothea_union  pearson   top_100                        0.7110                           0.2110                 10
   GSE126030_T_cells       blood        trrust_dorothea_union  pearson   top_500                        0.5934                           0.0934                 10
   GSE126030_T_cells       blood        trrust_dorothea_union spearman   top_100                        0.6160                           0.1160                 10
   GSE126030_T_cells       blood        trrust_dorothea_union spearman   top_500                        0.5494                           0.0494                 10
   GSE126030_T_cells bone_marrow                  dorothea_ab  pearson   top_100                        0.7410                           0.2410                 10
   GSE126030_T_cells bone_marrow                  dorothea_ab  pearson   top_500                        0.5862                           0.0862                 10
   GSE126030_T_cells bone_marrow                  dorothea_ab spearman   top_100                        0.7360                           0.2360                 10
   GSE126030_T_cells bone_marrow                  dorothea_ab spearman   top_500                        0.5998                           0.0998                 10
   GSE126030_T_cells bone_marrow                       trrust  pearson   top_100                        0.6450                           0.1450                 10
   GSE126030_T_cells bone_marrow                       trrust  pearson   top_500                        0.5614                           0.0614                 10
   GSE126030_T_cells bone_marrow                       trrust spearman   top_100                        0.6240                           0.1240                 10
   GSE126030_T_cells bone_marrow                       trrust spearman   top_500                        0.5440                           0.0440                 10
   GSE126030_T_cells bone_marrow trrust_dorothea_intersection  pearson   top_100                        0.6770                           0.1770                 10
   GSE126030_T_cells bone_marrow trrust_dorothea_intersection  pearson   top_500                        0.5240                           0.0240                 10
   GSE126030_T_cells bone_marrow trrust_dorothea_intersection spearman   top_100                        0.6380                           0.1380                 10
   GSE126030_T_cells bone_marrow trrust_dorothea_intersection spearman   top_500                        0.5196                           0.0196                 10
   GSE126030_T_cells bone_marrow        trrust_dorothea_union  pearson   top_100                        0.7440                           0.2440                 10
   GSE126030_T_cells bone_marrow        trrust_dorothea_union  pearson   top_500                        0.5962                           0.0962                 10
   GSE126030_T_cells bone_marrow        trrust_dorothea_union spearman   top_100                        0.7630                           0.2630                 10
   GSE126030_T_cells bone_marrow        trrust_dorothea_union spearman   top_500                        0.5888                           0.0888                 10
   GSE126030_T_cells        lung                  dorothea_ab  pearson   top_100                        0.7000                           0.2000                 10
   GSE126030_T_cells        lung                  dorothea_ab  pearson   top_500                        0.5824                           0.0824                 10
   GSE126030_T_cells        lung                  dorothea_ab spearman   top_100                        0.6930                           0.1930                 10
   GSE126030_T_cells        lung                  dorothea_ab spearman   top_500                        0.5758                           0.0758                 10
   GSE126030_T_cells        lung                       trrust  pearson   top_100                        0.6400                           0.1400                 10
   GSE126030_T_cells        lung                       trrust  pearson   top_500                        0.5302                           0.0302                 10
   GSE126030_T_cells        lung                       trrust spearman   top_100                        0.6590                           0.1590                 10
   GSE126030_T_cells        lung                       trrust spearman   top_500                        0.5186                           0.0186                 10
   GSE126030_T_cells        lung trrust_dorothea_intersection  pearson   top_100                        0.6200                           0.1200                 10
   GSE126030_T_cells        lung trrust_dorothea_intersection  pearson   top_500                        0.5142                           0.0142                 10
   GSE126030_T_cells        lung trrust_dorothea_intersection spearman   top_100                        0.5850                           0.0850                 10
   GSE126030_T_cells        lung trrust_dorothea_intersection spearman   top_500                        0.5192                           0.0192                 10
   GSE126030_T_cells        lung        trrust_dorothea_union  pearson   top_100                        0.7120                           0.2120                 10
   GSE126030_T_cells        lung        trrust_dorothea_union  pearson   top_500                        0.5796                           0.0796                 10
   GSE126030_T_cells        lung        trrust_dorothea_union spearman   top_100                        0.7080                           0.2080                 10
   GSE126030_T_cells        lung        trrust_dorothea_union spearman   top_500                        0.5734                           0.0734                 10
   GSE126030_T_cells  lymph_node                  dorothea_ab  pearson   top_100                        0.7240                           0.2240                 10
   GSE126030_T_cells  lymph_node                  dorothea_ab  pearson   top_500                        0.5902                           0.0902                 10
   GSE126030_T_cells  lymph_node                  dorothea_ab spearman   top_100                        0.7350                           0.2350                 10
   GSE126030_T_cells  lymph_node                  dorothea_ab spearman   top_500                        0.5818                           0.0818                 10
   GSE126030_T_cells  lymph_node                       trrust  pearson   top_100                        0.6730                           0.1730                 10
   GSE126030_T_cells  lymph_node                       trrust  pearson   top_500                        0.5304                           0.0304                 10
   GSE126030_T_cells  lymph_node                       trrust spearman   top_100                        0.6520                           0.1520                 10
   GSE126030_T_cells  lymph_node                       trrust spearman   top_500                        0.5338                           0.0338                 10
   GSE126030_T_cells  lymph_node trrust_dorothea_intersection  pearson   top_100                        0.6160                           0.1160                 10
   GSE126030_T_cells  lymph_node trrust_dorothea_intersection  pearson   top_500                        0.5180                           0.0180                 10
   GSE126030_T_cells  lymph_node trrust_dorothea_intersection spearman   top_100                        0.5940                           0.0940                 10
   GSE126030_T_cells  lymph_node trrust_dorothea_intersection spearman   top_500                        0.5074                           0.0074                  7
   GSE126030_T_cells  lymph_node        trrust_dorothea_union  pearson   top_100                        0.7520                           0.2520                 10
   GSE126030_T_cells  lymph_node        trrust_dorothea_union  pearson   top_500                        0.6076                           0.1076                 10
   GSE126030_T_cells  lymph_node        trrust_dorothea_union spearman   top_100                        0.7400                           0.2400                 10
   GSE126030_T_cells  lymph_node        trrust_dorothea_union spearman   top_500                        0.5780                           0.0780                 10
 GSE178429_IFN_gamma     ctrl_6h trrust_dorothea_intersection  pearson   top_100                        0.6420                           0.1420                 10
 GSE178429_IFN_gamma     ctrl_6h trrust_dorothea_intersection  pearson   top_500                        0.5152                           0.0152                 10
 GSE178429_IFN_gamma     ctrl_6h trrust_dorothea_intersection spearman   top_100                        0.6530                           0.1530                 10
 GSE178429_IFN_gamma     ctrl_6h trrust_dorothea_intersection spearman   top_500                        0.5122                           0.0122                 10
 GSE178429_IFN_gamma     ifng_6h trrust_dorothea_intersection  pearson   top_100                        0.5880                           0.0880                 10
 GSE178429_IFN_gamma     ifng_6h trrust_dorothea_intersection  pearson   top_500                        0.5502                           0.0502                 10
 GSE178429_IFN_gamma     ifng_6h trrust_dorothea_intersection spearman   top_100                        0.6230                           0.1230                 10
 GSE178429_IFN_gamma     ifng_6h trrust_dorothea_intersection spearman   top_500                        0.5652                           0.0652                 10
       Kang_IFN_beta        ctrl trrust_dorothea_intersection  pearson   top_100                        0.5380                           0.0380                 10
       Kang_IFN_beta        ctrl trrust_dorothea_intersection  pearson   top_500                        0.5052                           0.0052                 10
       Kang_IFN_beta        ctrl trrust_dorothea_intersection spearman   top_100                        0.5590                           0.0590                 10
       Kang_IFN_beta        ctrl trrust_dorothea_intersection spearman   top_500                        0.5084                           0.0084                  9
       Kang_IFN_beta        stim trrust_dorothea_intersection  pearson   top_100                        0.5470                           0.0470                 10
       Kang_IFN_beta        stim trrust_dorothea_intersection  pearson   top_500                        0.4980                          -0.0020                  3
       Kang_IFN_beta        stim trrust_dorothea_intersection spearman   top_100                        0.5230                           0.0230                  8
       Kang_IFN_beta        stim trrust_dorothea_intersection spearman   top_500                        0.5064                           0.0064                  9
NYGC_multimodal_PBMC         all                  dorothea_ab  pearson   top_100                        0.5780                           0.0780                 10
NYGC_multimodal_PBMC         all                  dorothea_ab  pearson   top_500                        0.5176                           0.0176                 10
NYGC_multimodal_PBMC         all                  dorothea_ab spearman   top_100                        0.5610                           0.0610                 10
NYGC_multimodal_PBMC         all                  dorothea_ab spearman   top_500                        0.5164                           0.0164                 10
NYGC_multimodal_PBMC         all                       trrust  pearson   top_100                        0.6350                           0.1350                 10
NYGC_multimodal_PBMC         all                       trrust  pearson   top_500                        0.5452                           0.0452                 10
NYGC_multimodal_PBMC         all                       trrust spearman   top_100                        0.6410                           0.1410                 10
NYGC_multimodal_PBMC         all                       trrust spearman   top_500                        0.5210                           0.0210                 10
NYGC_multimodal_PBMC         all trrust_dorothea_intersection  pearson   top_100                        0.5500                           0.0500                 10
NYGC_multimodal_PBMC         all trrust_dorothea_intersection  pearson   top_500                        0.5254                           0.0254                 10
NYGC_multimodal_PBMC         all trrust_dorothea_intersection spearman   top_100                        0.5510                           0.0510                 10
NYGC_multimodal_PBMC         all trrust_dorothea_intersection spearman   top_500                        0.5126                           0.0126                  9
NYGC_multimodal_PBMC         all        trrust_dorothea_union  pearson   top_100                        0.6170                           0.1170                 10
NYGC_multimodal_PBMC         all        trrust_dorothea_union  pearson   top_500                        0.5322                           0.0322                 10
NYGC_multimodal_PBMC         all        trrust_dorothea_union spearman   top_100                        0.5830                           0.0830                 10
NYGC_multimodal_PBMC         all        trrust_dorothea_union spearman   top_500                        0.5224                           0.0224                 10
             PBMC10k         all                  dorothea_ab  pearson   top_100                        0.5900                           0.0900                 10
             PBMC10k         all                  dorothea_ab  pearson   top_500                        0.5150                           0.0150                 10
             PBMC10k         all                  dorothea_ab spearman   top_100                        0.5900                           0.0900                 10
             PBMC10k         all                  dorothea_ab spearman   top_500                        0.5052                           0.0052                  9
             PBMC10k         all                       trrust  pearson   top_100                        0.6330                           0.1330                 10
             PBMC10k         all                       trrust  pearson   top_500                        0.5234                           0.0234                 10
             PBMC10k         all                       trrust spearman   top_100                        0.6490                           0.1490                 10
             PBMC10k         all                       trrust spearman   top_500                        0.5160                           0.0160                  9
             PBMC10k         all trrust_dorothea_intersection  pearson   top_100                        0.5540                           0.0540                 10
             PBMC10k         all trrust_dorothea_intersection  pearson   top_500                        0.5198                           0.0198                 10
             PBMC10k         all trrust_dorothea_intersection spearman   top_100                        0.5370                           0.0370                 10
             PBMC10k         all trrust_dorothea_intersection spearman   top_500                        0.5180                           0.0180                 10
             PBMC10k         all        trrust_dorothea_union  pearson   top_100                        0.6270                           0.1270                 10
             PBMC10k         all        trrust_dorothea_union  pearson   top_500                        0.5194                           0.0194                 10
             PBMC10k         all        trrust_dorothea_union spearman   top_100                        0.6280                           0.1280                 10
             PBMC10k         all        trrust_dorothea_union spearman   top_500                        0.5194                           0.0194                 10
