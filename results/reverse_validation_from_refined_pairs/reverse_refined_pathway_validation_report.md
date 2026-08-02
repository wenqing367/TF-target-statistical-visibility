# Reverse Pathway Validation From Refined Pairs

Scope: pathway validation is performed only on high-scoring pairs from the existing refined positive/background pair space.
The universe for each test is the target-gene set available in the same dataset, condition, edge set, and matching repeat.

- Repeat-level pathway tests: 42000
- Summary rows: 4200
- Analysis units with repeat universes: 28

Interpretation: pathway enrichment supports functional consistency of high-scoring targets; it does not prove causal regulation.

Top Pearson/Spearman pathway signals:
             dataset   condition              edge_set   metric top_level collection                         pathway  repeats_with_overlap  repeats_fdr_lt_0_05  mean_overlap_genes  mean_fold_enrichment  best_bh_fdr
   GSE126030_T_cells bone_marrow           dorothea_ab  pearson   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                19.9             11.896809 2.782440e-17
   GSE126030_T_cells bone_marrow trrust_dorothea_union spearman   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                17.8             11.220193 6.758272e-14
   GSE126030_T_cells       blood trrust_dorothea_union  pearson   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                16.8             11.182099 5.313170e-14
   GSE126030_T_cells bone_marrow trrust_dorothea_union  pearson   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                18.1             11.001631 1.765042e-15
   GSE126030_T_cells bone_marrow           dorothea_ab spearman   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                17.8             10.783289 8.094323e-15
   GSE126030_T_cells       blood           dorothea_ab  pearson   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                15.3             10.542853 5.136262e-12
   GSE126030_T_cells bone_marrow           dorothea_ab  pearson   top_100   Hallmark       Interferon Alpha Response                    10                   10                18.3             10.013026 4.912687e-15
   GSE126030_T_cells bone_marrow           dorothea_ab spearman   top_100   Hallmark       Interferon Alpha Response                    10                   10                17.3              9.587843 5.187763e-14
   GSE126030_T_cells bone_marrow trrust_dorothea_union spearman   top_100   Hallmark       Interferon Gamma Response                    10                   10                30.9              9.448338 7.927179e-22
   GSE126030_T_cells bone_marrow trrust_dorothea_union spearman   top_100   Hallmark       Interferon Alpha Response                    10                   10                17.0              9.394190 1.123230e-12
   GSE126030_T_cells bone_marrow trrust_dorothea_union  pearson   top_100   Hallmark       Interferon Alpha Response                    10                   10                17.5              9.321318 1.610249e-14
   GSE126030_T_cells bone_marrow trrust_dorothea_union  pearson   top_100   Hallmark       Interferon Gamma Response                    10                   10                30.0              8.850496 2.764828e-21
   GSE126030_T_cells bone_marrow           dorothea_ab  pearson   top_100   Hallmark       Interferon Gamma Response                    10                   10                27.5              8.000387 1.438606e-19
NYGC_multimodal_PBMC         all trrust_dorothea_union  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                12.0              7.753707 3.235060e-08
   GSE126030_T_cells bone_marrow           dorothea_ab  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 9.0              7.706254 7.407016e-06
   GSE126030_T_cells       blood trrust_dorothea_union spearman   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                12.0              7.672086 2.401494e-08
   GSE126030_T_cells       blood           dorothea_ab  pearson   top_100   Hallmark       Interferon Alpha Response                    10                   10                13.3              7.669090 2.132896e-09
   GSE126030_T_cells bone_marrow           dorothea_ab spearman   top_100   Hallmark       Interferon Gamma Response                    10                   10                25.7              7.574770 1.042412e-18
   GSE126030_T_cells bone_marrow trrust_dorothea_union spearman   top_100   Hallmark           Inflammatory Response                    10                   10                11.8              7.485297 1.891382e-07
   GSE126030_T_cells        lung trrust_dorothea_union  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 7.0              7.479537 1.068567e-04
   GSE126030_T_cells bone_marrow trrust_dorothea_union spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 8.0              7.421494 3.794557e-05
   GSE126030_T_cells       blood trrust_dorothea_union  pearson   top_100   Reactome            Interferon Signaling                    10                   10                24.9              7.320019 4.970589e-16
   GSE126030_T_cells       blood trrust_dorothea_union  pearson   top_100   Hallmark       Interferon Gamma Response                    10                   10                23.1              7.302469 6.467872e-14
   GSE126030_T_cells bone_marrow trrust_dorothea_union  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                11.8              7.223338 2.386498e-07
   GSE126030_T_cells bone_marrow trrust_dorothea_union  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 8.0              7.161599 4.672490e-05
   GSE126030_T_cells       blood           dorothea_ab spearman   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                10.9              7.159005 2.567127e-07
   GSE126030_T_cells  lymph_node trrust_dorothea_union  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                11.1              7.024190 3.802049e-07
   GSE126030_T_cells bone_marrow           dorothea_ab  pearson   top_100   Reactome            Interferon Signaling                    10                   10                24.7              6.979476 6.500788e-15
   GSE126030_T_cells        lung           dorothea_ab  pearson   top_100   Reactome Interferon Alpha Beta Signaling                    10                   10                10.3              6.825698 8.725000e-07
             PBMC10k         all trrust_dorothea_union spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 6.0              6.810055 1.693435e-05