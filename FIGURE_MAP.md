# Figure-to-code map

This map links EFP_analysis.R to the EFP paper figures. It records the manuscript/figure numbering reviewed during code preparation; the manuscript and figure deck are not runtime inputs. “Produces a panel” means the underlying plot is generated; manuscript-wide panel labels, cropping, manual labels and final layout in PowerPoint are separate steps.

This release retains the version 1.0.1 fix, which disables GAPIT's optional combined/circular GWAS plot bundle in
sections 08–09 to avoid the reproduced zero-step density-legend error. The
section 08 PC-comparison figures and section 10 rMVP paper-figure calls listed
below remain. Genotype diagnostics in section 07a and LD plots in sections
15–16 are unaffected by this plotting option.

**Numbering to confirm against the final publication:** the reviewed manuscript labels the PC comparison **S2** and the kinship heatmap **S3**; the slide deck labels kinship **S2** and PC comparison **S3**. Both are recorded below rather than silently selecting one numbering scheme.

| Paper panel | Code | Output / provenance | Limits or external additions |
|---|---|---|---|
| **1a** — genotype distributions | `efp_01()` | `outputs/Phenotypes/Xylose_barplot_ID_ordered_by_mean_faceted_by_Site.pdf/.png` | Mean within ID × Site × Breeder; independent within-site ordering. Extreme-genotype labels and min/max guide lines described in the caption are not implemented in the plotting functions. |
| **1b** — overlap/Venn diagram | **No code supplied** | External diagram | Archive the actual sample lists/counts and generating procedure. The release does not hard-code 344 or invent sample membership. |
| **1c** — site distributions | `efp_02()` | `Xylose_overlapping_histogram_by_Site1.pdf/.png` | 30 bins, explicit axis ticks. Uses mg/g xylose equivalents. The embedded histogram in the slide deck says mg/kg; reconcile the figure label with the reported units. |
| **1d** — site correlation | `efp_03()`; statistic from `efp_04()` | `Correlation_of_genotype_performance_between_Limagrain_RAGT.pdf/.png`, correlation CSV, `Heritability.csv` | The source plot/embedded image displays R², whereas the caption discusses Pearson r. The release labels both r and R² and exports P and N; heritability is exported separately rather than hard-coded on the plot. |
| **2a–c** — PCA pairs | GAPIT genotype diagnostics from `efp_07a()` and the archived GAPIT source | GAPIT-generated PCA graphics/tables | The study script delegates generation to GAPIT. It does not explicitly set the diagnostic PCA settings or name the exact output file chosen for assembly. Confirm and archive the actual source/output used. GWAS `PCA.total=0` controls covariates, not a requirement to show zero PCs in a diagnostic plot. |
| **2d–e** — power/type-I/FDR | `efp_07b()` | `outputs/GWAS/Compare_Power_FDR_test/`, GAPIT PDF and numeric CSV outputs | 100 replicates, h²=0.9, five QTNs retained. Source model vector includes SUPER (eight models); figure shows seven. Do not claim these are automatically the identical historical curves without the matching GAPIT source and simulation provenance. |
| **2f** — MLM inflation vs PCs | `efp_08()` | `PC_vs_genomic_inflation_only_MLM.pdf` | Full GWAS P values, lambda from median chi-square/median null chi-square. |
| **2g** — discoveries vs PCs, Average | `efp_08()` | `PC_vs_Bonferroni_significant_SNPs_only_Average.pdf` | BLINK/FarmCPU/MLMM; threshold 0.05/6791. The embedded slide plot uses the old 9.57e-6 label. Recompute/reconcile before final assembly. |
| **3a** — Manhattan | `efp_09()` → `efp_10()` | Per-model rMVP `Manhattan` JPG files and combined result CSVs | Assemble BLINK/FarmCPU/MLMM as desired. Manual peak labels/vertical guides are external. Release threshold is 5.1329637, not 5.02. |
| **3b** — QQ | `efp_10()` | Per-model rMVP `QQplot` JPG files | Final model/trait panel arrangement is assembled outside the code. |
| **4** — homozygous allele effects | `efp_11()` → `efp_12()` | `Boxplot_SigMarkers_alleles_vs_Xylose_homozygous_alleles_only_with_stats.pdf` | One marker per page, three trait facets. Choose seven lead-marker pages for the paper. ANOVA/Tukey at 0.05; means, percentage differences and pairwise adjusted P values are exported. Percentage/chromosome labels added in the final layout are not drawn automatically. |
| **5a** — EFP allele frequencies | `efp_13()` | `allele_frequencies_in_EFP_for_7_marker.pdf` | Requires externally curated percentage table 4. Bar values/markers should match the seven selected loci. Numerical labels in the assembled figure are external to the supplied plotting code. |
| **5b** — HFA count groups | `efp_14()` | `Boxplot_Average_Xylose_by_HFA_count_7_markers.pdf` and ANOVA/Tukey CSVs | Uses curated table 7; observed classes 0–6. Tukey letters use adjusted alpha 0.05 in the source, while the caption gives 0.001. Overall ANOVA P<0.001 is a separate statement. Percentage gain annotation is external. |
| **6** — lead-marker LD envelopes | `efp_16()` / `plot_marker_ld()` | `Fig6_LD_r2_0.3.pdf/.png`; comparison `Fig6_LD_r2_0.5.pdf/.png`; interval and linked-SNP CSVs | Seven figure markers; all requested markers get individual plots. Threshold is fixed per output; no automatic relaxation. Candidate-gene counts require external annotation/expression data and are not inserted. Intervals can change relative to older calculations; verify against final tables. |
| **7** — candidate-gene time course | `efp_17()` | `strongest_candidate_gene_from_each_QTL_expression_corrected.pdf/.png` plus summary CSV | Six selected genes, seven stages, mean ± SE. TimeCourse2 pooling follows source; n<2 has no estimated SE. |
| **S1** — marker QC | GAPIT diagnostics, especially `efp_07a()` | GAPIT-generated QC plots | Exact filename and contents depend on the archived GAPIT version. The pre-GWAS filter audit is also saved by 05/06 but is not the historical S1 graphic. |
| **Kinship: manuscript S3 / slides S2** | GAPIT diagnostics | GAPIT kinship heatmap | No independent heatmap command in the study script. Confirm the historical GAPIT output and kinship settings used. |
| **PC comparison: manuscript S2 / slides S3** | `efp_08()` | `PC_vs_genomic_inflation_lambda_without_MLM.pdf`, `PC_vs_number_of_Bonferroni_significant_SNPs_without_MLM.pdf` | Three multi-locus models across traits; assembled panels may differ in sizing. |
| **S4** — other germplasm frequencies | `efp_13()` | `allele_frequencies_by_7_marker.pdf` | EFP/Variety/Landrace values come from curated table 4. |
| **S5** — historical genome-wide LD | GAPIT diagnostic associated with genotype processing | Original GAPIT figure | The source does not provide a standalone historical S5 call/settings. Preserve that GAPIT output. `plot_genome_ld()` is explicitly a new adjacent-marker Pearson diagnostic and must not be relabelled as an exact historical reproduction. |
| **S6** — physical-position chromosome LD | `efp_15()` | `outputs/LD_full_chromosome/LD_heatmap_chr<CHR>.pdf` | This layout matches the physical-position heatmaps in the reviewed figure deck. Select chromosomes 1B, 3A, 5B, 6B, 7A for the paper. The triangular auxiliary images are additional outputs. |
| **S7** — original chromosome LD decay | `efp_15()` | `outputs/LD_full_chromosome/LD_decay_chr<CHR>.pdf`; bin CSVs | Starts at max(5 × median SNP spacing, 50 kb), targets median occupied-bin count ≥30, retains sparse bins. This differs from section 16; bin tables record what was used. |
| **Alternative S6/S7 outputs** | `efp_16()` | `outputs/GWAS/LD_simple/FigS6_LD_heatmaps.*`, `FigS7_LD_decay.*` | These names are retained from the section 16 implementation. Heatmaps are SNP-spaced; only decay bins individually containing ≥30 pairs are plotted. Distinguish them from section 15 in the final captions. |
| **S8** — expVIP interval-gene heatmaps in slide deck | **No generating code supplied** | External Wheat Expression Browser / expVIP output | Fig. 7 time-course code does not generate these log₂ TPM heatmaps. Record exact browser export, gene list, transformation and sample metadata. |

## Lead markers used for the seven Fig. 6 panels

`JB1B_662960818`, `AX-643813054`, `AX-94909517`, `AX-95679056`, `AX-643821508`, `AX-643812441`, `PER1`.

The source uses this explicit set. The significance table may contain additional SNPs at the same loci; the code does not silently replace the selected panel markers or automatically implement the manuscript's 1-Mb locus-grouping rule.
