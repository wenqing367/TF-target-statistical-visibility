# Reverse Pathway Validation From Refined Pairs

Scope: pathway validation is performed only on high-scoring pairs from the existing refined positive/background pair space.
The universe for each test is the target-gene set available in the same dataset, condition, edge set, and matching repeat.

- Repeat-level pathway tests: 18000
- Summary rows: 1800
- Analysis units with repeat universes: 12

Interpretation: pathway enrichment supports functional consistency of high-scoring targets; it does not prove causal regulation.

Top Pearson/Spearman pathway signals:
             dataset condition                     edge_set   metric top_level collection                         pathway  repeats_with_overlap  repeats_fdr_lt_0_05  mean_overlap_genes  mean_fold_enrichment  best_bh_fdr
NYGC_multimodal_PBMC       all        trrust_dorothea_union  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                12.0              7.753707 3.354528e-08
             PBMC10k       all        trrust_dorothea_union spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 6.0              6.810055 1.558908e-05
NYGC_multimodal_PBMC       all                  dorothea_ab  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                10.8              6.592177 6.246037e-07
             PBMC10k       all        trrust_dorothea_union  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 5.6              6.512222 1.051960e-04
             PBMC10k       all        trrust_dorothea_union  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                10.4              6.474984 5.587924e-07
             PBMC10k       all                       trrust spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 7.4              6.451583 4.325419e-05
NYGC_multimodal_PBMC       all        trrust_dorothea_union spearman   top_100   Hallmark           Inflammatory Response                    10                   10                10.4              6.389477 4.819925e-07
             PBMC10k       all                       trrust  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 7.0              6.167301 4.625763e-05
             PBMC10k       all        trrust_dorothea_union spearman   top_100   Hallmark           Inflammatory Response                    10                   10                10.0              6.078932 4.258756e-06
NYGC_multimodal_PBMC       all        trrust_dorothea_union spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 4.3              5.467904 3.395134e-03
NYGC_multimodal_PBMC       all                  dorothea_ab spearman   top_100   Hallmark           Inflammatory Response                    10                   10                 9.3              5.453904 8.641100e-06
             PBMC10k       all                  dorothea_ab  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                 8.6              5.190241 2.210938e-05
             PBMC10k       all                       trrust  pearson   top_100   Hallmark                      Complement                    10                   10                12.0              5.152576 5.814370e-08
             PBMC10k       all        trrust_dorothea_union  pearson   top_500   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                16.9              5.146047 1.349512e-09
             PBMC10k       all        trrust_dorothea_union spearman   top_500   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                16.8              5.099090 1.523172e-09
             PBMC10k       all                       trrust spearman   top_100   Hallmark                      Complement                    10                   10                11.8              5.011089 3.907512e-07
             PBMC10k       all                  dorothea_ab spearman   top_100   Hallmark           Inflammatory Response                    10                   10                 8.3              4.943354 1.414057e-04
NYGC_multimodal_PBMC       all                       trrust spearman   top_100   Hallmark           Inflammatory Response                    10                   10                11.5              4.908134 2.138282e-06
             PBMC10k       all        trrust_dorothea_union  pearson   top_100   Hallmark                      Complement                    10                   10                 8.8              4.895971 2.393003e-06
             PBMC10k       all                  dorothea_ab  pearson   top_500   Reactome Interferon Alpha Beta Signaling                    10                   10                18.1              4.813648 1.248709e-08
NYGC_multimodal_PBMC       all                       trrust  pearson   top_100   Hallmark           Inflammatory Response                    10                   10                11.5              4.775035 3.233511e-06
NYGC_multimodal_PBMC       all                       trrust  pearson   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 5.6              4.577925 2.668636e-03
             PBMC10k       all        trrust_dorothea_union  pearson   top_500   Reactome Interferon Alpha Beta Signaling                    10                   10                16.6              4.511232 8.422963e-08
             PBMC10k       all                       trrust  pearson   top_100   Hallmark       Interferon Gamma Response                    10                   10                16.5              4.413511 5.704056e-07
NYGC_multimodal_PBMC       all        trrust_dorothea_union  pearson   top_500   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                13.8              4.395779 1.707873e-06
             PBMC10k       all                  dorothea_ab  pearson   top_100   Hallmark                      Complement                    10                   10                 8.1              4.382541 1.323085e-05
NYGC_multimodal_PBMC       all                       trrust spearman   top_100   Hallmark        IL-6/JAK/STAT3 Signaling                    10                   10                 5.2              4.374171 2.775225e-03
             PBMC10k       all                  dorothea_ab spearman   top_500   Reactome Interferon Alpha Beta Signaling                    10                   10                16.3              4.351002 3.294803e-07
NYGC_multimodal_PBMC       all        trrust_dorothea_union  pearson   top_100   Hallmark       Interferon Gamma Response                    10                   10                10.4              4.346312 9.178843e-05
 GSE178429_IFN_gamma   ifng_6h trrust_dorothea_intersection spearman   top_100   Hallmark           Inflammatory Response                    10                   10                13.6              4.332035 5.767575e-08