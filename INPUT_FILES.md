# EFP input and output files

Paths below are relative to `project_root`. **No biological input data are bundled.** Schemas describe the file readers in EFP_analysis.R. The author reports a successful run of all sections; the biological input files were not attached for this packaging review. Preserve exact sample/marker names, including punctuation. CSVs have a header, no extra spreadsheet title rows, and no row-number column. Decimal separator is a dot. Use `NA` or an empty field for missing numeric observations; use the stated HapMap missing codes for genotypes. Do not use zero as a missing genotype or phenotype value.

## Required files

The first ten rows are externally prepared biological input files. The final row is bundled software. Five biological inputs retain their established paths under `outputs/`; do not delete them when clearing generated results. The previous `New_filtered_...` input is not required. The repository `.gitignore` explicitly allows these ten data files.

| File | Format and required fields | Used by |
|---|---|---|
| `inputs/EliteFibrePanel_phenotyping_for_plotting.csv` | Long CSV. `ID`: genotype ID; `Site`: `RAGT` or `Limagrain`; `Breeder`: `Limagrain`, `RAGT`, `F1 seed`/`F1 Seeds`, or `Saaten Union`; `Xylose`: numeric xylose equivalents, mg/g. One row per observation actually analysed. Record the level of aggregation and field/assay replication separately. | 01, 02, 04 |
| `inputs/EliteFibrePanel_phenotyping.csv` | Wide CSV, one unique `ID` per row. Numeric `RAGT`, `LIMAGRAIN`, `Average`; `Average` is the study's supplied across-site mean. It is not recomputed by the plotting functions. | 03, 12 |
| `inputs/EliteFibrePanel_phenotyping.txt` | Tab- or whitespace-delimited table with header. First column: unique sample ID (normally `Taxa` or `ID`). Remaining required columns: numeric `LIMAGRAIN`, `RAGT`, `Average`. IDs must match HapMap sample headers exactly. The runner orders phenotype rows by matched genotype IDs. | 08, 09 |
| `inputs/EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt` | Tab-delimited single-base HapMap; 11 metadata columns followed by sample columns. Calls: A/C/G/T, heterozygous IUPAC R/Y/S/W/K/M, missing N/NA/-/empty/dot. | 05 |
| `inputs/EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt` | Tab-delimited double-base HapMap; same metadata and sample layout. Homozygotes AA/CC/GG/TT; heterozygotes AC/CA/AG/GA/AT/TA/CG/GC/CT/TC/GT/TG; missing NN/NA/--/empty/dot. | 06 |
| `outputs/GWAS/AlleleFrequency/1_SigSNPs_to_calculate_allele_frequency.csv` | Manually prepared CSV: unique line `ID`, `group`, then marker columns. Calls are classified strings such as `HFA`, `LFA`, `Het`, `NN`, not nucleotide genotypes. The source does not show how high/low alleles were assigned or how external germplasm calls were matched. | 13 |
| `outputs/GWAS/AlleleFrequency/4_SigMarkers_allele_frequency_for_plotting.csv` | Curated long CSV: `Marker`, `GeneticResource` (`EFP`, `Variety`, `Landrace`), `Allele` (`HFA`, `LFA`), numeric `AlleleFrequency` in **percent, 0–100**. HFA/LFA pairs should sum to 100 within marker/resource after excluding heterozygous and missing calls, as stated in the Methods. | 13 |
| `outputs/GWAS/AlleleFrequency/5_for_HFA_counting_per_EFP_line.csv` | CSV: unique `ID`, followed only by the seven marker-call columns encoded as HFA/LFA/Het/NN or missing. Do not insert numeric phenotype columns here. | 14 |
| `outputs/GWAS/AlleleFrequency/7_HFA_counts_per_EFP_line_with_average_Xylose.csv` | Curated CSV with numeric `HFA_count` (integer 0–6 in the study) and `Average_Xylose` (mg/g). Retain `ID` for auditing the phenotype join. This file is not created by counting alone; the join to phenotypes occurred outside the supplied code. | 14 |
| `outputs/GWAS/Candidate_genes/Expression_time_course_data.csv` | CSV, one RNA-seq sample per row. Unique `sample`; `TimeCourse1`; `TimeCourse2` (study/series metadata); six gene columns of raw, nonnegative TPM. Exact gene IDs and stage labels below. | 17 |
| `vendor/gapit_functions.txt` | Included author-supplied GAPIT source, reporting version 4.1 (2026.5.19). Preserved byte-for-byte. Used with `source()`. An explicitly configured `EFP_GAPIT_SOURCE` path may override this location. | 07a, 07b, 08, 09 |

For all file paths above, “supplied” means provided by the repository author, even where the historical path begins with `outputs/`.

### HapMap details

Rows are markers; columns 12 onward are samples. The standard metadata sequence is `rs#`, `alleles`, `chrom`, `pos`, `strand`, `assembly#`, `center`, `protLSID`, `assayLSID`, `panelLSID`, `QCcode`. The first four fields provide marker ID, allele definition, chromosome and base-pair position. The reader preserves `rs#` and original sample names. At least two distinct homozygous classes must be present per retained marker. Heterozygosity and missingness thresholds are each 10%, with all sample columns in the denominator. The separate GAPIT MAF filter is 0.05.

The A/B/D subgenomes are represented as diploid, locus-specific markers in this analysis. Do not substitute hexaploid allele dosage 0–6 into the numerical genotype file.


**Note on the genotype input files**

To reduce file sizes for upload to GitHub, `EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt` and `EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt` include genotypic data for **chromosomes 1A and 1B only**. These are reduced demonstration datasets; the analyses reported in the paper used the complete genome-wide genotype datasets.

The reduced files illustrate the required input formats and can be used to test applicable parts of the workflow. They **cannot reproduce the full results reported in the paper**, and sections requiring markers on other chromosomes cannot be completed using these subsets. The significance threshold specified in the code relates to the full study dataset and has not been recalculated for these reduced files.

Reproducing the published analyses requires replacing these subsets with the complete genotype datasets provided in the publication’s supplementary tables.

### Expression stages and genes

`TimeCourse1` must be one of: `04 dpa`, `6-10 dpa`, `11-15 dpa`, `16-20 dpa`, `21-25 dpa`, `26-30 dpa`, `35 dpa`. Plot positions are respectively 4, 8, 13, 18, 23, 28 and 35 days post anthesis; these are display midpoints for grouped stages.

Required TPM columns:

- `TraesCS1B02G424500` — chromosome 1B, XTH/GH16.
- `TraesCS3A02G439700` — chromosome 3A, RING/U-box protein.
- `TraesCS5B02G356300` — chromosome 5B, UGPase.
- `TraesCS5B02G422600` — chromosome 5B, PPR-containing protein.
- `TraesCS6B02G042500` — chromosome 6B, PER1.
- `TraesCS7A02G349000` — chromosome 7A, TBL/DUF231 protein.

Means and SEs are calculated across rows within gene and stage; `TimeCourse2` is retained but not used to separate trajectories. Thus the interpretation assumes those rows are appropriate independent observations to pool. For n<2, SE is reported as missing and no confidence ribbon is implied. No RNA-seq alignment, transcript quantification, batch correction, or differential-expression analysis is performed by this script.

## Generated inputs consumed by later sections

| File or pattern | Produced by → read by | Content |
|---|---|---|
| `inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_sb_hmp_thresholds_0.1_0.1.txt` | 05 → 07a, 08, 09 | Filtered SB HapMap. May be supplied directly to skip filtering. |
| `inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_db_hmp_thresholds_0.1_0.1.txt` | 06 → 11 | Filtered DB HapMap. |
| `outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.Numerical.txt` | 07a → 07b, 15, 16 | Header; first column unique sample IDs; remaining columns SNP IDs; sample rows; calls 0/1/2/NA. LD rejects fractional dosages. |
| Same folder: `GAPIT.Genotype.map.txt` | 07a → 07b, 15, 16 | Header with `SNP`, `Chromosome`, `Position` (integer bp). IDs match numerical columns. LD accepts 1A–7D or 1–21 in order 1A,1B,1D,2A,…,7D. |
| `outputs/GWAS/results_<MODEL>_MAF_0.05_PCA_total_0/GAPIT.Association.GWAS_Results.<MODEL>.<TRAIT>(NYC).csv` | 09 → 10, 11 | Full results. `SNP`, `Chr`/`Chromosome`, `Pos`/`Position`, `P.value`/`P.Value`. Other GAPIT columns retained. Three models × three traits. Version-dependent filenames without `(NYC)` are also recognised. |
| `outputs/GWAS/PC_empirical_comparison/<MODEL>/PC_<0..10>/…GWAS_Results…csv` | 08 → 08 with `reuse_existing=TRUE` | Full results for MLM, MLMM, FarmCPU, BLINK and all three traits. Reuse only with unchanged inputs/software/settings. |
| `outputs/GWAS/list_of_sigMarkers.csv` | 11 → 16 | CSV with `sigMarker` column. All seven hard-coded figure markers must be present. Extra significant SNPs get separate plots. |
| `outputs/GWAS/transposed_hapmap_with_only_sigMarkers_data.csv` | 11 → 12 | **New consistent schema:** `ID` and DB genotype columns; one row per real sample, no metadata rows. Recreate using section 11 before running section 12 with older outputs. |
| `outputs/GWAS/LD_simple/LD_chr<CHR>.rds` | 16 → 16 with `resume=TRUE` | Saved LD map, r², complete-case counts, pair table and input calls. Older cache files are supported with an explicit message that exact allele-call verification is unavailable. |

## Output inventory by section

Every named output is relative to the same project root. GAPIT/rMVP also generate version-dependent filenames; preserve their complete output folders and logs.

| Section | Output directory | Files |
|---|---|---|
| 01 | `outputs/Phenotypes/` | `Xylose_barplot_ID_ordered_by_mean_faceted_by_Site.pdf/.png`; `Fig1a_genotype_site_means.csv`; `Fig1a_groups_without_phenotype.csv` |
| 02 | Same | `Xylose_overlapping_histogram_by_Site1.pdf/.png`; `Fig1c_rows_without_phenotype.csv` |
| 03 | Same | `Correlation_of_genotype_performance_between_Limagrain_RAGT.pdf/.png`; `Fig1d_correlation_statistics.csv` |
| 04 | Same | `Heritability.csv`, `Heritability_denominator_comparison.csv`, `Site_CV.csv`, `Heritability_model.rds` |
| 05–06 | `inputs/` | `Filtered_<raw stem>_thresholds_0.1_0.1.txt/.csv`; `Excluded_markers_filtered_<stem>_thresholds_0.1_0.1.txt`; `Filter_QC_<stem>_thresholds_0.1_0.1.csv` |
| 07a | `outputs/GWAS/Convert_HapMap_format_to_numerical/` | Numerical/map TXT files, `Logs_Convert_HapMap_format_to_numerical.txt`, GAPIT diagnostics |
| 07b | `outputs/GWAS/Compare_Power_FDR_test/` | `Logs_Power_compare.txt`, `Power_compare_returned_object.rds`, structure TXT, tabular return elements if available, GAPIT power/FDR/type-I CSV/PDF files |
| 08 | `outputs/GWAS/PC_empirical_comparison/` | `PC_empirical_comparison_summary.csv`; eight `PC_vs_*.pdf` plots; full GWAS result/log subfolders; `GAPIT_input_files.csv`, `GAPIT_phenotype_check.csv`, `GAPIT_chromosome_counts.csv`. A `models`/`pcs` subset uses `diagnostic_subset/` and does not create the full paper comparison plots. |
| 09 | `outputs/GWAS/results_<MODEL>_MAF_0.05_PCA_total_0/` | `Log_GWAS.txt`, full per-trait GWAS tables, GAPIT graphics/diagnostics |
| 10 | `outputs/GWAS/Combined_Manhattan_plots/<MODEL>/` | `Combined_GWAS_results_<MODEL>.csv`; rMVP JPG outputs identified by `Density`, `QQplot`, `Manhattan`, `Circular` memos |
| 11 | `outputs/GWAS/` | `sigMarkers_info.csv`, `list_of_sigMarkers.csv`, `transposed_hapmap_with_only_sigMarkers_data.csv` |
| 12 | `outputs/GWAS/SigMarkers_plots/` | `Boxplot_SigMarkers_alleles_vs_Xylose.pdf`; `_homozygous_alleles_only.pdf`; `_homozygous_alleles_only_with_stats.pdf`; `Allele_stats_ANOVA_meanDiff_homozygous_only.csv`; `Allele_Tukey_comparisons.csv`; merged wide/long phenotype/allele CSVs; unmatched-ID CSV |
| 13 | `outputs/GWAS/AlleleFrequency/` | `2_Allele_counts_per_group_per_SNP.csv`; `allele_frequencies_in_EFP_for_7_marker.pdf`; `allele_frequencies_by_7_marker.pdf` |
| 14 | Same | `6_HFA_counts_per_EFP_line.csv`; `Boxplot_Average_Xylose_by_HFA_count_7_markers.pdf`; `ANOVA_Tukey_letters_Average_Xylose_by_HFA_count.csv`; `HFA_Tukey_comparisons.csv` when estimable |
| 15 | `outputs/LD_full_chromosome/` | `LD_pairs_chr<CHR>.csv`; `LD_heatmap_chr<CHR>.pdf`; `LD_triangular_R2_chr<CHR>.pdf/.jpg`; `LD_decay_chr<CHR>.pdf`; `LD_decay_bins_chr<CHR>.csv` |
| 16 | `outputs/GWAS/LD_simple/` | `Fig6_LD_r2_0.3.pdf/.png`, `Fig6_LD_r2_0.5.pdf/.png`; `FigS6_LD_heatmaps.pdf/.png`, `FigS7_LD_decay.pdf/.png`; per-marker LD PDF/PNG and `_all_SNPs.csv`, `_linked_SNPs.csv`, `_interval.csv`; per-chromosome/local heatmaps; decay/count plots and bin CSVs; `LD_intervals_all_markers.csv`, `LD_pairs_chr<CHR>.csv.gz`, `LD_chr<CHR>.rds`, `Genotype_QC.csv`, `Sorted_map.csv`, `sessionInfo.txt` |
| 16 optional | Same | `FigS5_adjacent_diagnostics_REANALYSIS.pdf/.png`; `S5_adjacent_Pearson_pairs.csv` |
| 17 | `outputs/GWAS/Candidate_genes/Expression_time_course_plots11/` | `strongest_candidate_gene_from_each_QTL_expression_corrected.pdf/.png`; `strongest_candidate_gene_from_each_QTL_expression_corrected_summary.csv` |
| Runner | `outputs/Reproducibility/` | `sessionInfo.txt`, `installed_packages.csv`, `analysis_settings.R.txt`, `last_run_sections.csv`; GAPIT source checksum CSV when loaded |

## External preparation not implemented in EFP_analysis.R

EFP_analysis.R does not create the overlap/Venn diagram, extract annotation records for the LD intervals, calculate the gene counts overlaid in Fig. 6, classify HFA/LFA calls from raw nucleotides for external collections, perform the missing frequency/phenotype joins, download or quantify RNA-seq, generate the expVIP browser heatmaps, or perform automatic 1-Mb locus clustering. Archive the relevant prepared tables and describe these steps with their actual tools/versions and data access dates. Do not infer or hard-code new gene counts from a changed LD interval.
