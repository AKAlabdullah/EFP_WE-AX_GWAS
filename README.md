# EFP_WE-AX_GWAS
R code for **Genetic architecture of soluble arabinoxylan fibre in elite genotypes of bread wheat revealed by genome-wide association analysis**, by Abdul Kader Alabdullah and colleagues.
The analysis covers water-extractable arabinoxylan (WE-AX), expressed as xylose equivalents in wholemeal flour, in the Elite Fibre Panel: phenotypes, SNP filtering, GWAS, allele effects, linkage disequilibrium and candidate-gene expression.

**Release 1.0.3** contains one self-contained analysis script, one installer and the bundled GAPIT source. The earlier working R script is not included and is not required. The author reports successful execution of all sections. This packaging review preserves the executable analysis statements; it does not independently reproduce the study results.

## Run in R or RStudio

Open the extracted repository folder as your working directory. The two R files, this README, `inputs/`, `outputs/` and `vendor/` share the same project root. Add the ten study input files listed in [INPUT_FILES.md](INPUT_FILES.md) at their exact paths. Biological data are not included in this reviewed code archive.

```r
project_root <- "."

# Install once; omit these lines on later runs.
source("install_packages.R")
efp_install_packages(project_root)

# Load functions and generate the phenotype figures and summaries.
source("EFP_analysis.R")
efp_run(project_root, c("01", "02", "03", "04"))
```

Sourcing EFP_analysis.R only loads functions. `efp_run()` or an individual section call generates outputs. Output directories are created automatically; the required externally prepared tables must already be present. Rerunning a section replaces its named results. Sections 05–06 write filtered derivatives into `inputs/` without replacing the primary files.

```r
efp_run(project_root, "all")             # Entire workflow, in order
efp_run(project_root, "08")              # PC comparison only
efp_12(project_root)                     # Fig. 4 allele effects
results <- efp_16(project_root)          # Fig. 6 and LD outputs
efp_16(project_root, resume = TRUE)      # Replot unchanged saved LD
efp_17(project_root)                     # Fig. 7 expression
```

Section 07 has two IDs: **07a** (conversion/diagnostics) and **07b** (simulation). The `"all"` run includes both, all ten inputs, 100 simulation replicates and 44 model/PC combinations; it can take over 24 hours. Individual sections can run when their documented prerequisite files exist. Examples for every section appear at the end of EFP_analysis.R.

## Files in this repository

| File or folder | Purpose |
|---|---|
| `EFP_analysis.R` | Main script, settings, section functions and run examples |
| `install_packages.R` | One-time CRAN/Bioconductor dependency installation |
| `INPUT_FILES.md` | Required source files, generated intermediates, outputs and schemas |
| `FIGURE_MAP.md` | Figure panels, filenames and external assembly steps |
| `vendor/gapit_functions.txt` | Included GAPIT source, internally reporting version 4.1 |
| `vendor/gapit_README.md` and licence texts | GAPIT provenance, attribution and licence notices |
| `inputs/EliteFibrePanel_phenotyping_for_plotting.csv` | input file
| `inputs/EliteFibrePanel_phenotyping.csv` | Input file
| `inputs/EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt` | Input file
| `inputs/EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt` | Input file
| `outputs/GWAS/AlleleFrequency/1_SigSNPs_to_calculate_allele_frequency.csv` | Input file (manually prepared CSV) 
| `outputs/GWAS/AlleleFrequency/4_SigMarkers_allele_frequency_for_plotting.csv` | Input file (manually prepared CSV)
| `outputs/GWAS/AlleleFrequency/5_for_HFA_counting_per_EFP_line.csv` | Input file (manually prepared CSV)
| `outputs/GWAS/AlleleFrequency/7_HFA_counts_per_EFP_line_with_average_Xylose.csv` | Input file (manually prepared CSV)
| `outputs/GWAS/Candidate_genes/Expression_time_course_data.csv` | Input file (manually prepared CSV)

| Section | Analysis | Direct packages / software | Figure connection |
|---|---|---|---|
| 01 | Ranked genotype phenotypes | dplyr, ggplot2, tidytext | 1a |
| 02 | Site histograms | ggplot2 | 1c |
| 03 | Site correlation | ggplot2, scales; base stats | 1d |
| 04 | Heritability and site CV | lme4, dplyr | Statistic associated with 1d |
| 05 | Single-base HapMap QC | Base R | GWAS preprocessing |
| 06 | Double-base HapMap QC | Base R | Allele-effect preprocessing |
| 07a | HapMap conversion / genotype diagnostics | GAPIT source and dependencies | 2a–c; S1; kinship; S5, depending on archived GAPIT outputs |
| 07b | Model power, FDR and type I error | GAPIT source and dependencies | 2d–e |
| 08 | PCs 0–10 across four models | GAPIT, dplyr, ggplot2 | 2f–g and PC supplementary figure |
| 09 | Final GWAS | GAPIT source and dependencies | Results used in 3 |
| 10 | Manhattan, QQ, density and circular plots | rMVP | 3 |
| 11 | Significant SNPs and sample allele matrix | dplyr | Tables used in 4 and 6 |
| 12 | Allele effects / ANOVA / Tukey tests | dplyr, tidyr, ggplot2, gridExtra, multcompView | 4 |
| 13 | Favourable-allele counts and frequencies | dplyr, tidyr, ggplot2 | 5a, S4 |
| 14 | Favourable-allele count groups | dplyr, tidyr, readr, tibble, ggplot2, multcompView | 5b |
| 15 | Original chromosome heatmap / decay layout | snpStats, ggplot2, scales | S6–S7 in supplied slides |
| 16 | Lead-marker intervals and corrected LD workflow | snpStats, ggplot2, gridExtra, scales | 6; alternative S6–S7 |
| 17 | Expression time courses | dplyr, tidyr, readr, tibble, ggplot2 | 7 |


**Note on the genotype input files**

To reduce file sizes for upload to GitHub, `EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt` and `EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt` include genotypic data for **chromosomes 1A and 1B only**. These are reduced demonstration datasets; the analyses reported in the paper used the complete genome-wide genotype datasets.

The reduced files illustrate the required input formats and can be used to test applicable parts of the workflow. They **cannot reproduce the full results reported in the paper**, and sections requiring markers on other chromosomes cannot be completed using these subsets. The significance threshold specified in the code relates to the full study dataset and has not been recalculated for these reduced files.

Reproducing the published analyses requires replacing these subsets with the complete genotype datasets provided in the publication’s supplementary tables.

## Packages and GAPIT

Use R 4.3 or later and ggplot2 >= 3.4.0; record the actual environment from your successful run. `install_packages.R` lists direct analysis packages and the bundled GAPIT dependencies. snpStats and multtest are installed through Bioconductor. Base packages are supplied with R, and a standard R installation normally already provides the recommended package MASS. The installer only installs missing packages, except that ggplot2 is updated if below 3.4.0. It is not a historical package-version lockfile.

The bundled GAPIT file reports **2026.5.19, GAPIT 4.1**. It is preserved byte-for-byte; its checksum and upstream licence/citation details are in [vendor/gapit_README.md](vendor/gapit_README.md). No download is required when it is present. An existing `EFP_GAPIT_SOURCE` variable overrides the default path; use `Sys.unsetenv("EFP_GAPIT_SOURCE")` to clear a previous override. The loader records which file it used.

## Section 08 and reproducibility records

Section 08 and the final GWAS in section 09 both use the filtered single-base HapMap produced by section 05. A separate `New_filtered_...` input is not required. The PC comparison writes to `outputs/GWAS/PC_empirical_comparison/`. Its default `reuse_existing=FALSE` reruns GWAS. Only use `reuse_existing=TRUE` for complete results known to match the same inputs, settings and GAPIT source: this attached implementation checks result-file presence, not cache provenance.

Sections 08–09 retain `Multiple_analysis=FALSE` to avoid the optional GAPIT circular-plot error. Section 08 still produces PC-comparison PDFs, and section 10 produces the rMVP plots. For a diagnostic subset:

```r
efp_08(project_root, models = "MLM", pcs = 0)
```

This writes under `diagnostic_subset/` and does not replace the full comparison plots. GAPIT logs capture warnings, error calls and session details. Input checks save file checksums, phenotype summaries and chromosome counts.

`efp_run()` records `sessionInfo.txt`, `installed_packages.csv`, `analysis_settings.R.txt` and `last_run_sections.csv` in `outputs/Reproducibility/`. These files describe the **latest runner invocation** and are overwritten by later calls; `last_run_sections.csv` is not a cumulative record of separately run sections. Preserve the successful run records before further runs.

## Statistical settings and figure scope

The default threshold is `0.05 / 6791 = 7.3626859e-6` (−log10 P = 5.1329637). This is based on the total marker count, not an estimated effective count. Final GWAS uses MAF 0.05 and zero PCs. The attached release retains these settings.

Heritability uses `Xylose ~ Site + (1 | ID)` and exports both analysed-observation and historical row-count denominator summaries. It does not fit genotype-by-site variance. LD sections 15 and 16 use different plot layouts and decay bins; select the matching figure captions. The optional section 16 adjacent-marker diagnostic is a new analysis and does not reproduce historical Fig. S5. See [FIGURE_MAP.md](FIGURE_MAP.md) for manuscript checks and external figure steps.

## References

- [GAPIT project and version-specific citations](https://github.com/jiabowang/GAPIT).
- [snpStats](https://bioconductor.org/packages/snpStats/) and [rMVP](https://cran.r-project.org/package=rMVP).
