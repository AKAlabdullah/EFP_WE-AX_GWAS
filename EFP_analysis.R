# EFP analysis code for GitHub
# Genetic architecture of soluble arabinoxylan fibre in elite genotypes of
# bread wheat revealed by genome-wide association analysis
# Author of the study code: Abdul Kader Alabdullah, John Innes Centre
# Self-contained publication script; no earlier EFP script is required.
# Release 1.0.3 (9 October 2026): repository documentation and packaging review.
#
# RUN: source("EFP_analysis.R"); efp_run("path/to/project", c("01", "02", "03", "04"))
# Each efp_XX() function is independent when its documented files are available.
# Sourcing this file LOADS functions; efp_run() or efp_XX() GENERATES outputs.
# Examples for every section appear at the end. No RStudio API is required.
# See INPUT_FILES.md for exact schemas and FIGURE_MAP.md for figure provenance.
# Methods and result-sensitive corrections are recorded in RELEASE_NOTES.md.

# Study settings for this release; see RELEASE_NOTES.md for interpretation.
efp_settings <- list(
  maf = 0.05, pc_total = 0L, n_markers = 6791L, alpha = 0.05,
  models = c("MLMM", "FarmCPU", "BLINK"),
  traits = c("LIMAGRAIN", "RAGT", "Average"),
  power_models = c("GLM", "MLM", "MLMM", "CMLM", "ECMLM", "FarmCPU", "BLINK", "SUPER"),
  power_nrep = 100L, power_h2 = 0.9, power_nqtn = 5L,
  # Reproducibility seed used by the simulation, final GWAS and plot jitter.
  seed = 20261008L,
  hapmap = "inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_sb_hmp_thresholds_0.1_0.1.txt",
  # PC comparison and final GWAS use the same filtered single-base HapMap.
  pc_hapmap = "inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_sb_hmp_thresholds_0.1_0.1.txt"
)

# Small shared helpers. All analyses below remain ordinary R functions.
.efp_root <- function(path) {
  if (!dir.exists(path)) stop("Project directory not found: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
.efp_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages using install_packages.R: ", paste(missing, collapse = ", "))
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  if ("ggplot2" %in% pkgs && utils::packageVersion("ggplot2") < "3.4.0")
    stop("ggplot2 >= 3.4.0 is required.")
}
.efp_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create output folder: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
.efp_csv <- function(path, columns = character()) {
  if (!file.exists(path)) stop("Missing input file: ", path, " (see INPUT_FILES.md)")
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN"))
  if (anyDuplicated(names(x))) stop("Duplicate column names in ", path)
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(path, " requires columns: ", paste(absent, collapse = ", "))
  x
}
.efp_number <- function(x, label) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.na(x) & !is.finite(y))) stop("Non-numeric/non-finite values in ", label)
  y
}
.efp_unique_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
    stop(label, ": IDs must be nonempty and unique.")
}
.efp_hapmap <- function(path, header = TRUE) {
  if (!file.exists(path)) stop("Missing HapMap: ", path)
  # comment.char='' is essential: the standard first header is rs#.
  h <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = c("", "NA", "NaN"))
  if (ncol(h) < 13L || nrow(h) == 0) stop("HapMap requires 11 metadata columns and >=2 samples: ", path)
  .efp_unique_ids(h[[1]], paste(path, "SNP"))
  .efp_unique_ids(names(h)[-(1:11)], paste(path, "sample"))
  if (header) return(h)
  # GAPIT expects the header row inside the G object (header=FALSE convention).
  ans <- rbind(names(h), as.matrix(h)); colnames(ans) <- NULL
  as.data.frame(ans, stringsAsFactors = FALSE)
}
.efp_gwas_inputs <- function(hapmap, report_dir = NULL) {
  G <- .efp_hapmap(hapmap, header = FALSE)
  # Validate the HapMap passed by the PC comparison or final GWAS section.
  position <- suppressWarnings(as.numeric(G[-1, 4]))
  bad <- !is.finite(position) | position < 0 | position != floor(position)
  if (any(bad)) stop("Invalid base-pair positions in ", hapmap, "; SNPs: ",
    paste(head(G[-1, 1][bad], 10), collapse = ", "))
  chromosome <- as.character(G[-1, 3])
  wheat_chr <- as.vector(t(outer(1:7, c("A", "B", "D"), paste0)))
  if (anyNA(chromosome) || any(!chromosome %in% c(wheat_chr, as.character(1:21))))
    stop("HapMap chromosome column must contain 1A..7D or numeric 1..21: ", hapmap)
  # Validate without changing allele coding, marker order, or chromosome codes.
  calls <- as.matrix(G[-1, -(1:11), drop = FALSE])
  valid <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "N",
    "AA", "CC", "GG", "TT", "AC", "CA", "AG", "GA", "AT", "TA", "CG", "GC",
    "CT", "TC", "GT", "TG", "NN", "-", "--", ".")
  invalid <- setdiff(unique(as.vector(calls[!is.na(calls)])), valid)
  if (length(invalid)) stop("Unexpected HapMap calls in ", hapmap, ": ",
    paste(head(invalid, 10), collapse = ", "))
  path <- "inputs/EliteFibrePanel_phenotyping.txt"
  if (!file.exists(path)) stop("Missing phenotype table: ", path)
  Y <- read.table(path, header = TRUE, check.names = FALSE, comment.char = "",
                  stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA", "NaN", "."))
  if (!all(efp_settings$traits %in% names(Y))) stop("GWAS phenotype table requires LIMAGRAIN, RAGT, Average.")
  .efp_unique_ids(as.character(Y[[1]]), "GWAS phenotype")
  for (tr in efp_settings$traits) Y[[tr]] <- .efp_number(Y[[tr]], tr)
  Y <- Y[, c(names(Y)[1], efp_settings$traits), drop = FALSE]
  ids <- as.character(unlist(G[1, -(1:11)]))
  if (length(setdiff(Y[[1]], ids))) stop("Phenotype IDs absent from HapMap. Match original IDs; do not replace punctuation.")
  Y <- Y[order(match(Y[[1]], ids)), , drop = FALSE]
  qc <- do.call(rbind, lapply(efp_settings$traits, function(tr) {
    y <- Y[[tr]]; finite <- is.finite(y)
    data.frame(Trait = tr, N_observed = sum(finite), N_missing = sum(!finite),
      N_distinct = length(unique(y[finite])),
      Minimum = if (any(finite)) min(y[finite]) else NA_real_,
      Maximum = if (any(finite)) max(y[finite]) else NA_real_)
  }))
  if (!is.null(report_dir)) {
    report_dir <- .efp_dir(report_dir)
    write.csv(qc, file.path(report_dir, "GAPIT_phenotype_check.csv"), row.names = FALSE)
    write.csv(data.frame(File = c(normalizePath(hapmap, winslash = "/"),
      normalizePath(path, winslash = "/")), MD5 = unname(tools::md5sum(c(hapmap, path))),
      Rows = c(nrow(G) - 1L, nrow(Y)), Columns = c(ncol(G), ncol(Y))),
      file.path(report_dir, "GAPIT_input_files.csv"), row.names = FALSE)
    write.csv(as.data.frame(table(Chromosome = chromosome)),
      file.path(report_dir, "GAPIT_chromosome_counts.csv"), row.names = FALSE)
  }
  if (any(qc$N_observed < 3L | qc$N_distinct < 2L))
    stop("Each GWAS trait needs >=3 observed values and phenotype variation. Check the phenotype table.")
  list(Y = Y, G = G)
}
.efp_gapit <- function(project_root) {
  path <- Sys.getenv("EFP_GAPIT_SOURCE", unset = file.path(project_root, "vendor", "gapit_functions.txt"))
  if (!file.exists(path)) stop("Supply the study's archived GAPIT source at vendor/gapit_functions.txt. See README.md.")
  # GAPIT is sourced software in this project; it need not be an installed R package.
  source(path, local = .GlobalEnv, encoding = "UTF-8")
  if (!exists("GAPIT", mode = "function")) stop("GAPIT() not defined by ", path)
  provenance <- .efp_dir(file.path(project_root, "outputs", "Reproducibility"))
  version <- if (exists("GAPIT.0000", mode = "function"))
    tryCatch(as.character(GAPIT.0000()), error = function(e) "Not reported") else "Not reported"
  write.csv(data.frame(File = normalizePath(path, winslash = "/"),
            MD5 = unname(tools::md5sum(path)), GAPIT_version = paste(version, collapse = "; ")),
            file.path(provenance, "GAPIT_source_used.csv"), row.names = FALSE)
  invisible(path)
}
.efp_log <- function(path, code) {
  # Record errors BEFORE stack unwinding. The outer runner otherwise hides
  # the internal GAPIT call. Warnings are logged immediately, not only the first 50.
  # Restore sinks, directory, warning option and devices on success or failure.
  force(path); old_dir <- getwd(); old_sinks <- sink.number(); old_message <- sink.number(type = "message")
  old_devices <- grDevices::dev.list(); old_active <- grDevices::dev.cur()
  old_warn <- getOption("warn"); failed <- FALSE
  log_path <- file.path(normalizePath(dirname(path), winslash = "/"), basename(path))
  con <- file(path, open = "wt")
  sink(con); sink(con, type = "message")
  on.exit({
    sink(type = "message")
    if (old_message != 2L) sink(getConnection(old_message), type = "message")
    while (sink.number() > old_sinks) sink()
    close(con); setwd(old_dir); options(warn = old_warn)
    for (d in setdiff(grDevices::dev.list(), old_devices)) grDevices::dev.off(d)
    if (old_active > 1 && old_active %in% grDevices::dev.list()) grDevices::dev.set(old_active)
    if (failed) message("GAPIT diagnostic log: ", log_path)
  }, add = TRUE)
  cat("Started: ", format(Sys.time()), "\nR: ", R.version.string, "\n\n", sep = "")
  withCallingHandlers(force(code),
    warning = function(w) {
      cat("\nWARNING: ", conditionMessage(w), "\nCall: ",
        paste(deparse(conditionCall(w)), collapse = " "), "\n", sep = "")
      # Do not suppress warnings or change GAPIT's warning/error policy.
    },
    error = function(e) {
      failed <<- TRUE
      calls <- sys.calls(); frames <- sys.frames()
      cat("\n===== GAPIT ERROR =====\n", conditionMessage(e), "\nCall: ",
        paste(deparse(conditionCall(e)), collapse = " "), "\n\nCall stack:\n", sep = "")
      for (i in seq_along(calls))
        cat(i, ": ", paste(deparse(calls[[i]]), collapse = " "), "\n", sep = "")
      # For the reported seq() error, record scalar sequence arguments only.
      # Do not dump whole genotype matrices or replace base R's seq().
      for (i in seq_along(calls)) {
        if (grepl("^seq(\\.default)?\\(", paste(deparse(calls[[i]]), collapse = " "))) {
          cat("\nSequence frame ", i, ":\n", sep = "")
          for (nm in c("from", "to", "by", "length.out", "n")) {
            value <- tryCatch(get(nm, envir = frames[[i]], inherits = FALSE), error = function(e) NULL)
            if (is.numeric(value) && length(value) == 1L)
              cat(nm, " = ", format(value, digits = 17), "\n", sep = "")
          }
        }
      }
      cat("\nSession:\n"); print(sessionInfo())
      flush(con)
    })
}
.efp_result_file <- function(run_dir, model, trait) {
  # File naming differs between GAPIT versions. Prefer exact NYC results.
  candidates <- file.path(run_dir, paste0("GAPIT.Association.GWAS_Results.", model, ".", trait,
                                          c("(NYC).csv", ".csv")))
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[1])
  files <- list.files(run_dir, pattern = "\\.csv$", full.names = TRUE)
  found <- files[grepl(paste0(model, ".", trait), basename(files), fixed = TRUE) &
                   grepl("GWAS[._]Results", basename(files)) & !grepl("Kansas", basename(files))]
  if (length(found) != 1) stop("Need one full GWAS result for ", model, "/", trait, " in ", run_dir)
  found
}
.efp_results <- function(file) {
  x <- .efp_csv(file)
  # Harmonise aliases without losing or recalculating GAPIT's other columns.
  aliases <- list(SNP = c("SNP"), Chr = c("Chr", "Chromosome"), Pos = c("Pos", "Position"),
                  P.value = c("P.value", "P.Value", "p.value", "P", "p", "pvalue"))
  for (target in names(aliases)) {
    hit <- intersect(aliases[[target]], names(x))
    if (!length(hit)) stop("Missing ", target, " in ", file)
    names(x)[match(hit[1], names(x))] <- target
  }
  x$SNP <- as.character(x$SNP); .efp_unique_ids(x$SNP, file)
  x$Pos <- .efp_number(x$Pos, "GWAS position")
  x$P.value <- .efp_number(x$P.value, "P value")
  if (any(x$P.value < 0 | x$P.value > 1, na.rm = TRUE)) stop("P values outside [0,1].")
  x
}


# Stage exports before replacing the final file, and reject empty outputs.
.efp_ggsave <- function(filename = NULL, plot, file = NULL, ...) {
  force(plot)
  if (is.null(filename)) filename <- file
  if (is.null(filename)) stop("A plot filename is required.")
  folder <- .efp_dir(dirname(filename)); ext <- tools::file_ext(filename)
  tmp <- tempfile(tmpdir=folder,fileext=paste0(".",ext))
  on.exit(unlink(tmp),add=TRUE)
  ggplot2::ggsave(tmp,plot=plot,device=ext,...)
  if (!file.exists(tmp) || file.info(tmp)$size < 100) stop("Empty plot: ",filename)
  if (!file.copy(tmp,filename,overwrite=TRUE) || file.info(filename)$size != file.info(tmp)$size)
    stop("Could not save plot: ",filename)
  invisible(filename)
}
.efp_multipage <- function(plots, filename, width=4, height=4) {
  if (!length(plots)) stop("No plots available for ",filename)
  previous <- grDevices::dev.cur()
  grDevices::pdf(NULL); build_device <- grDevices::dev.cur()
  on.exit({
    if (build_device %in% grDevices::dev.list()) grDevices::dev.off(build_device)
    if (previous > 1 && previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  },add=TRUE)
  # Build before opening the output PDF, avoiding a blank first page.
  grobs <- lapply(plots,ggplot2::ggplotGrob)
  grDevices::dev.off(build_device)
  tmp <- tempfile(tmpdir=.efp_dir(dirname(filename)),fileext=".pdf")
  on.exit(unlink(tmp),add=TRUE)
  grDevices::pdf(tmp,width=width,height=height,onefile=TRUE)
  output_device <- grDevices::dev.cur()
  on.exit(if(output_device %in% grDevices::dev.list()) grDevices::dev.off(output_device),add=TRUE)
  for (g in grobs) { grid::grid.newpage(); grid::grid.draw(g) }
  grDevices::dev.off(output_device)
  if (!file.exists(tmp) || file.info(tmp)$size < 100) stop("Empty PDF: ",filename)
  if (!file.copy(tmp,filename,overwrite=TRUE) || file.info(filename)$size != file.info(tmp)$size)
    stop("Could not save PDF: ",filename)
  invisible(filename)
}

# =============================================================================
# 01 | Ranked phenotypes
# Packages: dplyr, ggplot2, tidytext.
# Inputs: inputs/EliteFibrePanel_phenotyping_for_plotting.csv; ID, Site, Breeder, Xylose (mg/g).
# Outputs: outputs/Phenotypes/Xylose_barplot_ID_ordered_by_mean_faceted_by_Site.{pdf,png}.
# Tables: Fig1a_genotype_site_means.csv; Fig1a_groups_without_phenotype.csv.
# Paper: Fig. 1a; manual highlights and assembly are external.
# Export reduced from 50 x 20 to 20 x 8 inches to avoid very large PNGs.
# Analysis section: 01 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_01 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "ggplot2", "tidytext"))
  output_dir1 <- .efp_dir("outputs/Phenotypes")

  XyloseData <- .efp_csv("inputs/EliteFibrePanel_phenotyping_for_plotting.csv", c("ID", "Site", "Breeder", "Xylose"))
  XyloseData$Breeder[XyloseData$Breeder == "F1 seed"] <- "F1 Seeds"

  trait <- "Xylose"

  # Ensure numeric
  XyloseData[[trait]] <- .efp_number(XyloseData[[trait]], trait)

  # Summarise to mean per ID per Site
  Xylose_summary <- XyloseData %>%
    group_by(ID, Site, Breeder) %>%
    summarise(
      mean_Xylose = mean(Xylose, na.rm = TRUE),
      .groups = "drop"
    )

  # A genotype/site group with no measurements has mean NaN; it cannot form a bar.
  missing_groups <- Xylose_summary[!is.finite(Xylose_summary$mean_Xylose), ]
  write.csv(missing_groups, file.path(output_dir1, "Fig1a_groups_without_phenotype.csv"), row.names = FALSE)
  if (nrow(missing_groups)) message("Section 01: ", nrow(missing_groups),
    " groups have no phenotype; listed in Fig1a_groups_without_phenotype.csv.")
  if (!any(is.finite(Xylose_summary$mean_Xylose))) stop("No observed phenotypes for Fig. 1a.")

  p1 <- ggplot(
    Xylose_summary[is.finite(Xylose_summary$mean_Xylose), ],
    aes(
      x = reorder_within(ID, mean_Xylose, Site),
      y = mean_Xylose,
      fill = Breeder
    )
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.5),
      #colour = "black",
      width = 0.7
    ) +
    facet_wrap(~ Site, nrow = 2, scales = "free_x") +
    scale_x_reordered() +
    ggtitle("Xylose content by genotype, ordered from low to high") +
    xlab("Genotype ID") +
    ylab("Xylose equivalents (mg/g)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size = 10,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14),
      axis.line = element_line(color = "black"),
        axis.ticks = element_line(colour = "black"),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    scale_fill_manual(
      values = c(
        "Limagrain"      = "#c89600",
        "RAGT"           = "#00affa",
        "F1 Seeds"        = "#6a3d9a",
        "Saaten Union"   = "#33a02c"
      )
    )

  # Save the plot
  .efp_ggsave(
    file = paste0(output_dir1, "/", trait, "_barplot_ID_ordered_by_mean_faceted_by_Site.pdf"),
    plot = p1,
    dpi = 300,
    limitsize = FALSE,
    height = 8,
    width = 20
  )

  .efp_ggsave(
    file = paste0(output_dir1, "/", trait, "_barplot_ID_ordered_by_mean_faceted_by_Site.png"),
    plot = p1,
    dpi = 300,
    limitsize = FALSE,
    height = 8,
    width = 20
  )
  write.csv(Xylose_summary, file.path(output_dir1, "Fig1a_genotype_site_means.csv"), row.names = FALSE)
  invisible(Xylose_summary)
}


# =============================================================================
# 02 | Site histograms
# Packages: ggplot2.
# Inputs: same long phenotype CSV as 01; 30 bins, original observations.
# Outputs: outputs/Phenotypes/Xylose_overlapping_histogram_by_Site1.{pdf,png}.
# Missing phenotypes: Fig1c_rows_without_phenotype.csv (no imputation).
# Paper: Fig. 1c; units are mg/g (not mg/kg).
# Analysis section: 02 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_02 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("ggplot2"))
  output_dir1 <- .efp_dir("outputs/Phenotypes")

  XyloseData <- .efp_csv("inputs/EliteFibrePanel_phenotyping_for_plotting.csv", c("ID", "Site", "Xylose"))
 
  trait <- c("Xylose")
  # Convert traits to numeric if needed
  XyloseData[[trait]] <- .efp_number(XyloseData[[trait]], trait)
  missing_rows <- !is.finite(XyloseData[[trait]])
  write.csv(XyloseData[missing_rows, ], file.path(output_dir1, "Fig1c_rows_without_phenotype.csv"), row.names = FALSE)
  if (any(missing_rows)) message("Section 02: ", sum(missing_rows),
    " missing phenotypes excluded; listed in Fig1c_rows_without_phenotype.csv.")
  XyloseData <- XyloseData[!missing_rows, ]
  if (!nrow(XyloseData)) stop("No observed phenotypes for Fig. 1c.")
  p1 <- ggplot(XyloseData, aes(x = .data[[trait]], fill = Site, color = Site)) + 
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +  # Removed alpha for solid color
    ggtitle(paste("Histogram of", trait)) + 
    xlab("Xylose equivalents mg/g") + 
    ylab("Count") +
    theme_minimal() + 
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.text = element_text(size = 14),
      axis.line = element_line(color = "black"),
        axis.ticks = element_line(colour = "black"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.margin = unit(c(1, 3, 1, 1), "cm")
    ) +
    scale_fill_manual(values = c("Limagrain" = "#c89600", "RAGT" = "#00affa")) +
    scale_color_manual(values = c("Limagrain" = "#c89600", "RAGT" = "#00affa"))

  # Save the plot
  .efp_ggsave(file = paste0(output_dir1, "/", trait, "_overlapping_histogram_by_Site1.pdf"),
         plot = p1, dpi = 300, limitsize = FALSE, height = 3, width = 6)
  .efp_ggsave(file = paste0(output_dir1, "/", trait, "_overlapping_histogram_by_Site1.png"),
         plot = p1, dpi = 300, limitsize = FALSE, height = 3, width = 6)
}


# =============================================================================
# 03 | Between-site correlation
# Packages: ggplot2, scales.
# Inputs: inputs/EliteFibrePanel_phenotyping.csv; one ID per row, RAGT and LIMAGRAIN numeric.
# Outputs: outputs/Phenotypes/Correlation_of_genotype_performance_between_Limagrain_RAGT.{pdf,png}; Fig1d_correlation_statistics.csv.
# Paper: Fig. 1d; both Pearson r and its square are explicitly labelled.
# Analysis section: 03 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_03 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("ggplot2", "scales"))
  output_dir1 <- .efp_dir("outputs/Phenotypes")
  Mydata <- .efp_csv("inputs/EliteFibrePanel_phenotyping.csv", c("ID", "RAGT", "LIMAGRAIN"))
  .efp_unique_ids(Mydata$ID, "Wide phenotype table")
  Mydata$RAGT <- .efp_number(Mydata$RAGT, "RAGT")
  Mydata$LIMAGRAIN <- .efp_number(Mydata$LIMAGRAIN, "LIMAGRAIN")
  Mydata <- Mydata[complete.cases(Mydata[, c("RAGT", "LIMAGRAIN")]), ]
  if (nrow(Mydata) < 3 || sd(Mydata$RAGT) == 0 || sd(Mydata$LIMAGRAIN) == 0) stop("Need >=3 variable paired observations.")
  ct <- cor.test(Mydata$RAGT, Mydata$LIMAGRAIN, method = "pearson")
  fit <- lm(LIMAGRAIN ~ RAGT, data = Mydata)
  r <- unname(ct$estimate)
  label <- sprintf("r = %.2f; R² = %.2f", r, r^2)
  p2 <- ggplot(Mydata, aes(RAGT, LIMAGRAIN)) +
    geom_point(color = "#07873a", alpha = 0.4, size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red2") +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.05, vjust = 1.3, size = 4) +
    labs(x = "RAGT (xylose equivalents, mg/g)", y = "Limagrain (xylose equivalents, mg/g)") +
    scale_x_continuous(breaks = scales::breaks_pretty(5)) +
    scale_y_continuous(breaks = scales::breaks_pretty(5)) + theme_classic(base_size = 12)
  for (ext in c("pdf", "png")) .efp_ggsave(file.path(output_dir1, paste0("Correlation_of_genotype_performance_between_Limagrain_RAGT.", ext)), p2, width = 4, height = 4, dpi = 300)
  result <- data.frame(N_pairs = nrow(Mydata), Pearson_r = r, R_squared = r^2, P_value = ct$p.value,
                       Slope = unname(coef(fit)[2]), Intercept = unname(coef(fit)[1]))
  write.csv(result, file.path(output_dir1, "Fig1d_correlation_statistics.csv"), row.names = FALSE)
  invisible(result)
}


# =============================================================================
# 04 | Heritability and coefficients of variation
# Packages: lme4, dplyr.
# Inputs: long phenotype CSV from 01; complete ID/Site/Xylose observations.
# Outputs: outputs/Phenotypes/Heritability.csv, Site_CV.csv, Heritability_model.rds.
# Audit: Heritability_denominator_comparison.csv (analysed vs all input rows).
# Paper: Statistic associated with Fig. 1d; no figure generated.
# Analysis section: 04 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_04 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("lme4", "dplyr"))
  # Retains Xylose ~ Site + (1|ID), REML; Site is fixed.
  # H2 = Vg/(Vg + Ve/r), r = mean observations per ID (not field blocks).
  # This model does NOT separate genotype-by-site variance or estimate a GxE test.
  XyloseData <- .efp_csv("inputs/EliteFibrePanel_phenotyping_for_plotting.csv", c("ID", "Site", "Xylose"))
  XyloseData$Xylose <- .efp_number(XyloseData$Xylose, "Xylose")
  # Historical denominator included rows subsequently omitted by lmer().
  # Retain it only as an explicitly labelled comparison, not field replication.
  r_input_rows <- mean(table(XyloseData$ID))
  input_rows <- nrow(XyloseData)
  XyloseData <- XyloseData[complete.cases(XyloseData[, c("ID", "Site", "Xylose")]), ]

  # Fit mixed model: ID is random, Site is fixed
  model <- lmer(Xylose ~ Site + (1 | ID), data = XyloseData, REML = TRUE)

  # Extract variance components
  var_comp <- as.data.frame(VarCorr(model))
  gen_var <- var_comp$vcov[var_comp$grp == "ID"]
  res_var <- var_comp$vcov[var_comp$grp == "Residual"]

  # Estimate H² (broad-sense heritability)
  r <- mean(table(XyloseData$ID))  # average analysed observations per genotype
  H2 <- gen_var / (gen_var + (res_var / r))

  cat("Broad-sense heritability (H²):", round(H2, 3), "\n\n")

  cv_results <- XyloseData %>% group_by(Site) %>% summarise(N = n(), Mean = mean(Xylose), SD = sd(Xylose), CV = SD / Mean * 100, .groups = "drop")
  folder <- .efp_dir("outputs/Phenotypes")
  write.csv(data.frame(Genetic_variance = gen_var, Residual_variance = res_var, Mean_observations_per_ID = r, H2 = H2, Singular = lme4::isSingular(model)), file.path(folder, "Heritability.csv"), row.names = FALSE)
  write.csv(data.frame(Denominator = c("Analysed observations per ID", "All input rows per ID (historical comparison)"),
    Mean_observations_per_ID = c(r, r_input_rows),
    H2 = gen_var / (gen_var + res_var / c(r, r_input_rows)),
    Input_rows = input_rows, Analysed_rows = nrow(XyloseData),
    Omitted_rows = input_rows - nrow(XyloseData)),
    file.path(folder, "Heritability_denominator_comparison.csv"), row.names = FALSE)
  write.csv(cv_results, file.path(folder, "Site_CV.csv"), row.names = FALSE)
  saveRDS(model, file.path(folder, "Heritability_model.rds"))
  print(cv_results); invisible(list(H2 = H2, CV = cv_results))
}


# Shared filtering for sections 05 and 06; rates divide by all samples per SNP.
.efp_filter_hapmap <- function(path, encoding, het_thresh = 0.1, miss_thresh = 0.1) {
  h <- .efp_hapmap(path)
  g <- toupper(trimws(as.matrix(h[, -(1:11), drop = FALSE])))
  dim(g) <- c(nrow(h), ncol(h) - 11L)
  missing_codes <- if (encoding == "sb") c("N", "NA", "-", "", ".") else c("NN", "NA", "--", "", ".")
  g[is.na(g) | g %in% missing_codes] <- if (encoding == "sb") "N" else "NN"
  hom <- if (encoding == "sb") c("A", "C", "G", "T") else c("AA", "CC", "GG", "TT")
  het <- if (encoding == "sb") c("R", "Y", "S", "W", "K", "M") else c("AC", "CA", "AG", "GA", "AT", "TA", "CG", "GC", "CT", "TC", "GT", "TG")
  missing <- if (encoding == "sb") "N" else "NN"
  invalid <- setdiff(unique(as.vector(g)), c(hom, het, missing))
  if (length(invalid)) stop("Unexpected HapMap call(s): ", paste(invalid, collapse = ", "))
  logs <- lapply(seq_len(nrow(g)), function(i) {
    x <- g[i, ]; nhom <- length(unique(x[x %in% hom]))
    reasons <- c(if (length(unique(x[x != missing])) <= 1) "Monomorphic",
                 if (nhom < 2) "Only one homozygous allele",
                 if (mean(x %in% het) > het_thresh) "Heterozygosity exceeds threshold",
                 if (mean(x == missing) > miss_thresh) "Missingness exceeds threshold")
    data.frame(SNP = h[i, 1], Missing_rate = mean(x == missing), Heterozygosity = mean(x %in% het),
               Retained = !length(reasons), Reason = paste(reasons, collapse = " & "))
  })
  qc <- do.call(rbind, logs); h[, -(1:11)] <- g
  stem <- paste0(tools::file_path_sans_ext(basename(path)), "_thresholds_", het_thresh, "_", miss_thresh)
  write.table(h[qc$Retained, , drop = FALSE], paste0("inputs/Filtered_", stem, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(h[qc$Retained, , drop = FALSE], paste0("inputs/Filtered_", stem, ".csv"), row.names = FALSE)
  write.table(qc[!qc$Retained, ], paste0("inputs/Excluded_markers_filtered_", stem, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(qc, paste0("inputs/Filter_QC_", stem, ".csv"), row.names = FALSE)
  message(sum(qc$Retained), " of ", nrow(qc), " SNPs retained (MAF filter is applied separately by GAPIT).")
  invisible(qc)
}


# =============================================================================
# 05 | Single-base HapMap filtering
# Packages: base R only.
# Inputs: inputs/EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt; tab-delimited HapMap with 11 metadata columns.
# Outputs: inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_sb_hmp_thresholds_0.1_0.1.{txt,csv}; Excluded_markers_filtered_*.txt and Filter_QC_*.csv.
# Paper: Preprocessing for GWAS/allelic plots; not the GAPIT S1 figure.
# Analysis section: 05 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_05 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(character())
  .efp_filter_hapmap("inputs/EFP_genotypes_AxiomTaNG1_KASP_sb_hmp.txt", "sb")
}


# =============================================================================
# 06 | Double-base HapMap filtering
# Packages: base R only.
# Inputs: inputs/EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt; tab-delimited HapMap with 11 metadata columns.
# Outputs: inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_db_hmp_thresholds_0.1_0.1.{txt,csv}; Excluded_markers_filtered_*.txt and Filter_QC_*.csv.
# Paper: Preprocessing for GWAS/allelic plots; not the GAPIT S1 figure.
# Analysis section: 06 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_06 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(character())
  .efp_filter_hapmap("inputs/EFP_genotypes_AxiomTaNG1_KASP_db_hmp.txt", "db")
}


# =============================================================================
# 07a | Convert HapMap and generate GAPIT diagnostics
# Packages: GAPIT source plus its dependencies; see install_packages.R.
# Inputs: filtered SB HapMap from 05; vendor/gapit_functions.txt.
# Outputs: outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.{Numerical,map}.txt plus GAPIT diagnostics and log.
# Paper: GAPIT PCA/QC/kinship/LD diagnostics contribute to Figs. 2a-c, S1, kinship S2/S3 and S5; see FIGURE_MAP.md.
# Analysis section: 07a (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_07a <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(character())
  .efp_gapit(project_root)
  myG <- .efp_hapmap(efp_settings$hapmap, header = FALSE)
  folder <- .efp_dir("outputs/GWAS/Convert_HapMap_format_to_numerical")
  setwd(folder)
  fit <- .efp_log("Logs_Convert_HapMap_format_to_numerical.txt", GAPIT(G = myG, output.numerical = TRUE))
  # GAPIT's own outputs retain their filenames; require them before downstream use.
  needed <- c("GAPIT.Genotype.Numerical.txt", "GAPIT.Genotype.map.txt")
  if (any(!file.exists(needed))) stop("This GAPIT source did not write the expected numerical/map files: ", paste(needed, collapse = ", "))
  invisible(fit)
}


# =============================================================================
# 07b | Simulation of model power, FDR and type I error
# Packages: GAPIT source plus its dependencies; see install_packages.R.
# Inputs: numerical genotype/map files from 07a; GAPIT source.
# Outputs: outputs/GWAS/Compare_Power_FDR_test/; GAPIT PDFs/CSVs, returned-object RDS and log.
# Paper: Fig. 2d-e; supplied code runs eight models including SUPER, figure legend has seven.
# Analysis section: 07b (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_07b <- function(project_root = ".", nrep = efp_settings$power_nrep) {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(character())
  .efp_gapit(project_root)
  # 100 replicates can take >24 hours. nrep=2 is only a smoke test.
  myGD <- read.table(
    "./outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.Numerical.txt",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  myGM <- read.table(
    "./outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.map.txt",
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # Convert wheat chromosomes (1A–7D) into numbers (1–21)
  chr.old <- c("1A","1B","1D","2A","2B","2D","3A","3B","3D",
               "4A","4B","4D","5A","5B","5D","6A","6B","6D",
               "7A","7B","7D")
  chr.new <- as.character(1:21)

  # -------------------- Align GD and GM by SNP name (CRITICAL) --------------------
  # Taxa column is the FIRST column in GD (do not assume its name)
  taxa_col <- colnames(myGD)[1]

  # marker IDs
  gd_snps <- colnames(myGD)[-1]
  gm_snps <- as.character(myGM$SNP)

  # intersection (common markers)
  common <- intersect(gd_snps, gm_snps)

  # reorder BOTH to the same order
  myGD2 <- myGD[, c(taxa_col, common), drop = FALSE]
  colnames(myGD2)[1] <- "Taxa"   # GAPIT expects this name in many functions

  myGM2 <- myGM[match(common, gm_snps), , drop = FALSE]

  stopifnot(all(colnames(myGD2)[-1] == myGM2$SNP))  # MUST be TRUE now

  # -------------------- Chromosome filter/recode on aligned GM2 --------------------
  chr_text <- toupper(sub("^CHR", "", as.character(myGM2$Chromosome)))
  chr.idx <- match(chr_text, chr.old)
  num_idx <- suppressWarnings(as.integer(chr_text))
  chr.idx[is.na(chr.idx) & num_idx %in% 1:21] <- num_idx[is.na(chr.idx) & num_idx %in% 1:21]
  keep_chr <- !is.na(chr.idx)

  # keep only core chromosomes and recode
  myGM2 <- myGM2[keep_chr, ]
  myGM2$Chromosome <- as.numeric(chr.new[chr.idx[keep_chr]])

  # apply same keep_chr to aligned GD2 (drop the same SNPs)
  myGD2 <- cbind(
    myGD2[, 1, drop = FALSE],
    myGD2[, -1, drop = FALSE][, keep_chr, drop = FALSE]
  )

  # -------------------- MAF filter on aligned GD2 --------------------
  p   <- colMeans(myGD2[, -1, drop = FALSE], na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  keep_maf <- is.finite(maf) & maf > efp_settings$maf

  # Filter GD and GM consistently
  GD_filt <- cbind(
    myGD2[, 1, drop = FALSE],
    myGD2[, -1, drop = FALSE][, keep_maf, drop = FALSE]
  )
  GM_filt <- myGM2[keep_maf, ]

  # -------------------- Sanity checks --------------------
  stopifnot(ncol(GD_filt) - 1 == nrow(GM_filt))
  stopifnot(all(colnames(GD_filt)[-1] == GM_filt$SNP))

  cat("Individuals:", nrow(GD_filt), "\n")
  cat("Markers kept:", ncol(GD_filt) - 1, "\n")
  print(table(GM_filt$Chromosome))

  if (!nrow(GM_filt)) stop("No markers remain for simulation.")
  folder <- .efp_dir("outputs/GWAS/Compare_Power_FDR_test")
  setwd(folder)
  set.seed(efp_settings$seed)
  power.out <- .efp_log("Logs_Power_compare.txt", GAPIT.Power.compare(
   GD = GD_filt, GM = GM_filt, nrep = nrep, h2 = efp_settings$power_h2,
   model = efp_settings$power_models, NQTN = efp_settings$power_nqtn))
  saveRDS(power.out, "Power_compare_returned_object.rds")
  writeLines(capture.output(str(power.out)), "Power_compare_returned_object_structure.txt")
  # GAPIT also writes *.Power.by.FDR.*.csv / *.Power.by.TypeI.*.csv (version-dependent).
  # Retain these numeric tables alongside its PDF, including for reviewer responses.
  if (is.list(power.out)) for (nm in names(power.out)) {
   x <- power.out[[nm]]
   if (is.matrix(x) || is.data.frame(x)) write.csv(x, paste0("Power_compare_", gsub("[^A-Za-z0-9]", "_", nm), ".csv"), row.names = FALSE)
  }
  invisible(power.out)
}


# =============================================================================
# 08 | Empirical PC comparison
# Packages: dplyr, ggplot2, GAPIT source plus its dependencies.
# Inputs: GWAS phenotype TXT; filtered SB HapMap from 05; vendor GAPIT source.
# Outputs: outputs/GWAS/PC_empirical_comparison/; full runs, PC_empirical_comparison_summary.csv and PC_vs_*.pdf.
# Also: GAPIT input-file/phenotype/chromosome checks; detailed per-run error logs.
# Optional models/pcs subsets write under diagnostic_subset/ without paper plots.
# Paper: Fig. 2f-g and PC supplementary figure (S2 in manuscript / S3 in slide deck).
# Analysis section: 08 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_08 <- function(project_root = ".", reuse_existing = FALSE,
                   models = c("MLM", "BLINK", "FarmCPU", "MLMM"), pcs = 0:10) {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "ggplot2"))
  .efp_gapit(project_root)
  # reuse_existing=TRUE only with unchanged genotype, phenotype, GAPIT and settings.
  # --------------------------
  # Inputs
  # --------------------------
  # --------------------------
  # Settings
  # --------------------------
  allowed_models <- c("MLM", "BLINK", "FarmCPU", "MLMM")
  if (!length(models) || anyNA(models) || any(!models %in% allowed_models))
    stop("models must be selected from MLM, BLINK, FarmCPU, MLMM.")
  if (!is.numeric(pcs) || !length(pcs) || anyNA(pcs) || any(!pcs %in% 0:10))
    stop("pcs must be selected from the integers 0:10.")
  pcs_to_try <- unique(pcs)
  models_to_test <- unique(models)
  subset_run <- !setequal(models, allowed_models) || !setequal(pcs, 0:10)
  maf_threshold <- efp_settings$maf
  bonf_p <- efp_settings$alpha / efp_settings$n_markers

  out_dir <- "./outputs/GWAS/PC_empirical_comparison"
  # A diagnostic subset must not overwrite the full paper comparison/figures.
  if (subset_run) out_dir <- file.path(out_dir, "diagnostic_subset")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- .efp_gwas_inputs(efp_settings$pc_hapmap, report_dir = out_dir)
  myY <- inputs$Y; myG <- inputs$G
  message("PC-comparison HapMap: ", efp_settings$pc_hapmap)

  # --------------------------
  # Functions
  # --------------------------
  calc_lambda <- function(p) {
    p <- as.numeric(p)
    p <- p[!is.na(p) & p > 0 & p <= 1]
    if (length(p) < 10) return(NA_real_)
    chisq_obs <- qchisq(p, df = 1, lower.tail = FALSE)
    median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
  }

  gapit_output_exists <- function(run_dir) {
    all(vapply(efp_settings$traits, function(tr)
      tryCatch({ .efp_result_file(run_dir, model_name, tr); TRUE }, error = function(e) FALSE), logical(1)))
  }

  # --------------------------
  # Run GAPIT for each model and PC number
  # --------------------------
  summary_list <- list()
  counter <- 1

  for (model_name in models_to_test) {
  
    for (k in pcs_to_try) {
    
      cat("\n============================================================\n")
      cat("Model:", model_name, "| PCA.total =", k, "\n")
      cat("============================================================\n")
    
      run_dir <- file.path(out_dir, model_name, paste0("PC_", k))
      run_dir <- .efp_dir(run_dir)
    
      # ----------------------------------------------------------
      # Skip GAPIT if output files already exist
      # ----------------------------------------------------------
      if (reuse_existing && gapit_output_exists(run_dir)) {
      
        cat("Existing GAPIT output found. Skipping GWAS run:\n")
        cat(run_dir, "\n")
      
      } else {
      
        cat("Running GAPIT; existing results are reused only with reuse_existing=TRUE.\n")
      
        old_wd <- getwd()
        setwd(run_dir)
      
        log_file <- paste0("Log_", model_name, "_PC_", k, ".txt")
        # GAPIT's optional combined/circular plots can fail in Densitplot() at
        # seq(0, maxbin.num, len) when len=0, after the GWAS tables are written.
        # PC figures below use those tables; combined paper GWAS plots use section 10.
        # Keep result-file output and the statistical settings, omit only that plot bundle.
        fit <- .efp_log(log_file, GAPIT(Y = myY, G = myG, PCA.total = k,
          model = model_name, SNP.MAF = maf_threshold, Multiple_analysis = FALSE, file.output = TRUE))
        setwd(old_wd)
      }
    
      # --------------------------
      # Read full GAPIT result files for lambda and SNP counts
      # --------------------------
      full <- lapply(efp_settings$traits, function(tr) {
        x <- .efp_results(.efp_result_file(run_dir, model_name, tr))
        pv <- x$P.value; good <- is.finite(pv)
        data.frame(Model = model_name, PCA.total = k, Trait = tr,
          N_SNPs = sum(good), Lambda = calc_lambda(pv),
          Min_P = if (any(good)) min(pv[good]) else NA_real_,
          SNPs_p_Bonf = sum(pv <= bonf_p, na.rm = TRUE), Status = "OK")
      })
      summary_k <- do.call(rbind, full)
      summary_list[[counter]] <- summary_k
      counter <- counter + 1
    }
  }

  # --------------------------
  # Combine, filter and save summary
  # --------------------------
  pc_summary <- bind_rows(summary_list) %>%
    filter(!grepl("\\(Kansas\\)", Trait)) %>%
    mutate(
      Trait = trimws(gsub("\\s*\\(NYC\\)", "", Trait))
    )

  write.csv(
    pc_summary,
    file.path(out_dir, "PC_empirical_comparison_summary.csv"),
    row.names = FALSE
  )

  print(pc_summary)

  cat("\nFinished.\n")
  cat("Summary table saved to:\n")
  cat(file.path(out_dir, "PC_empirical_comparison_summary.csv"), "\n")
  if (subset_run) {
    message("Diagnostic subset complete. Run efp_08(project_root) for all models/PCs and paper plots.")
    return(invisible(pc_summary))
  }

  # ==========================================================================
  # Plot all models all traits 
  # ==========================================================================
  # --------------------------
  # Plot genomic inflation factor
  # --------------------------
  p_lambda0 <- ggplot(
    pc_summary %>% filter(Status == "OK"),
    aes(x = PCA.total, y = Lambda, group = Model, colour = Model)
  ) +
    geom_line() +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    labs(
      x = "Number of principal components included in GWAS",
      y = expression("Genomic inflation factor (" * lambda * ")"),
      title = "Effect of PCA correction on genomic inflation across GWAS models"
    ) +
    theme_minimal(base_size = 14)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_genomic_inflation_lambda.pdf"),
    p_lambda0,
    width = 12,
    height = 6
  )

  # --------------------------
  # Plot significant SNPs at Bonferroni threshold
  # --------------------------
  p_sig_bonf0 <- ggplot(
    pc_summary %>% filter(Status == "OK"),
    aes(x = PCA.total, y = SNPs_p_Bonf, group = Model, colour = Model)
  ) +
    geom_line() +
    geom_point(size = 2.5) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    labs(
      x = "Number of principal components included in GWAS",
      y = paste0("Number of SNPs with P <= ", format(bonf_p, digits = 3)),
      title = "Effect of PCA correction on Bonferroni-significant associations"
    ) +
    theme_minimal(base_size = 14)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_number_of_Bonferroni_significant_SNPs.pdf"),
    p_sig_bonf0,
    width = 12,
    height = 6
  )

  # ==========================================================================
  # Re-plot models without MLM all traits 
  # ==========================================================================

  ## Remove MLM and keep only Average
  pc_summary_withoutMLM <- pc_summary %>%
    filter(Model != "MLM")

  # --------------------------
  # Plot genomic inflation factor
  # --------------------------
  p_lambda1 <- ggplot(
    pc_summary_withoutMLM %>% filter(Status == "OK"),
    aes(x = PCA.total, y = Lambda, group = Model, colour = Model)
  ) +
    geom_line() +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    labs(
      x = "Number of principal components included in GWAS",
      y = expression("Genomic inflation factor (" * lambda * ")"),
      title = "Effect of PCA correction on genomic inflation across GWAS models"
    ) +
    theme_minimal(base_size = 14)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_genomic_inflation_lambda_without_MLM.pdf"),
    p_lambda1,
    width = 12,
    height = 4
  )

  # --------------------------
  # Plot significant SNPs at Bonferroni threshold
  # --------------------------
  p_sig_bonf1 <- ggplot(
    pc_summary_withoutMLM %>% filter(Status == "OK"),
    aes(x = PCA.total, y = SNPs_p_Bonf, group = Model, colour = Model)
  ) +
    geom_line() +
    geom_point(size = 2.5) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    labs(
      x = "Number of principal components included in GWAS",
      y = paste0("Number of SNPs with P <= ", format(bonf_p, digits = 3)),
      title = "Effect of PCA correction on Bonferroni-significant associations"
    ) +
    theme_minimal(base_size = 14)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_number_of_Bonferroni_significant_SNPs_without_MLM.pdf"),
    p_sig_bonf1,
    width = 12,
    height = 4
  )

  # ==========================================================================
  # Re-plot only model MLM (Fig. 2f)
  # ==========================================================================

  ## Keep only MLM
  pc_summary_MLM <- pc_summary %>%
    filter(grepl("MLM", Model))%>%
    filter(!grepl("MLMM", Model))

  trait_cols <- c(
    "LIMAGRAIN" = "#00A9E0",
    "RAGT"      = "#FDBB00",
    "Average"   = "#7030A0"
  )

  # ----------------------------------------------------
  # Plot genomic inflation
  # ----------------------------------------------------
  p_lambda2 <- ggplot(
    pc_summary_MLM %>% filter(Status == "OK"),
    aes(
      x = PCA.total,
      y = Lambda,
      colour = Trait,
      group = Trait
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_colour_manual(values = trait_cols) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    labs(
      x = "Number of principal components included in GWAS",
      y = expression("Genomic inflation factor (" * lambda * ")"),
      title = "Effect of PCA correction on genomic inflation",
      colour = "Trait"
    ) +
    theme_bw(base_size = 12)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_genomic_inflation_only_MLM.pdf"),
    p_lambda2,
    width = 6,
    height = 3
  )

  # ----------------------------------------------------
  # Plot Bonferroni significant SNPs
  # ----------------------------------------------------
  p_bonf2 <- ggplot(
    pc_summary_MLM %>% filter(Status=="OK"),
    aes(PCA.total,
        SNPs_p_Bonf,
        colour=Trait,
        group=Trait)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_colour_manual(values = trait_cols) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    labs(
      x="Number of principal components included in GWAS",
      y=paste0("Number of SNPs with P <= ", format(bonf_p, digits = 3)),
      title="Effect of PCA correction on Bonferroni-significant associations"
    ) +
    theme_bw(base_size=12)

  .efp_ggsave(
    file.path(out_dir,
              "PC_vs_Bonferroni_significant_SNPs_only_MLM.pdf"),
    p_bonf2,
    width=6,
    height=3
  )

  # ==========================================================================
  # Re-plot models without MLM only Average 
  # ==========================================================================
  ## Remove MLM and keep only Average
  pc_summary_withoutMLM_onlyAverage <- pc_summary_withoutMLM %>%
    filter(grepl("Average", Trait))

  # ----------------------------------------------------
  # Plot genomic inflation
  # ----------------------------------------------------
  p_lambda3 <- ggplot(
    pc_summary_withoutMLM_onlyAverage %>% filter(Status == "OK"),
    aes(
      x = PCA.total,
      y = Lambda,
      colour = Model,
      group = Model
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    labs(
      x = "Number of principal components included in GWAS",
      y = expression("Genomic inflation factor (" * lambda * ")"),
      title = "Effect of PCA correction on genomic inflation",
      colour = "Trait"
    ) +
    theme_bw(base_size = 12)

  .efp_ggsave(
    file.path(out_dir, "PC_vs_genomic_inflation_only_Average.pdf"),
    p_lambda3,
    width = 6,
    height = 3
  )

  # -----------------------------------------------------------------------------
  # Plot Bonferroni significant SNPs all models with only Average (Fig. 2g)
  # -----------------------------------------------------------------------------
  p_bonf3 <- ggplot(
    pc_summary_withoutMLM_onlyAverage %>% filter(Status=="OK"),
    aes(PCA.total,
        SNPs_p_Bonf,
        colour=Model,
        group=Model)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    labs(
      x="Number of principal components included in GWAS",
      y=paste0("Number of SNPs with P <= ", format(bonf_p, digits = 3)),
      title="Effect of PCA correction on Bonferroni-significant associations"
    ) +
    theme_bw(base_size=12)

  .efp_ggsave(
    file.path(out_dir,
              "PC_vs_Bonferroni_significant_SNPs_only_Average.pdf"),
    p_bonf3,
    width=6,
    height=3
  )
}


# =============================================================================
# 09 | Final GWAS
# Packages: GAPIT source plus its dependencies; see install_packages.R.
# Inputs: GWAS phenotype TXT; filtered SB HapMap; local GAPIT source.
# Outputs: outputs/GWAS/results_<model>_MAF_0.05_PCA_total_0/; full GAPIT results, diagnostics, Log_GWAS.txt.
# Paper: Source results for Fig. 3; PCA.total=0; BLINK, FarmCPU, MLMM.
# Analysis section: 09 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_09 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(character())
  .efp_gapit(project_root)
  x <- .efp_gwas_inputs(efp_settings$hapmap)
  for (model_name in efp_settings$models) {
   folder <- .efp_dir(file.path(project_root, "outputs/GWAS", paste0("results_", model_name,
     "_MAF_", efp_settings$maf, "_PCA_total_", efp_settings$pc_total)))
   setwd(folder); message("GWAS: ", model_name)
   set.seed(efp_settings$seed)
   # Section 10 produces the paper GWAS plots; avoid the optional GAPIT
   # combined/circular plotting failure documented in section 08.
   fit <- .efp_log("Log_GWAS.txt", GAPIT(Y = x$Y, G = x$G,
      PCA.total = efp_settings$pc_total, SNP.MAF = efp_settings$maf,
      Multiple_analysis = FALSE, model = model_name, Random.model = FALSE,
      file.output = TRUE))
  }
  invisible(NULL)
}


# =============================================================================
# 10 | GWAS plots for each model
# Packages: rMVP.
# Inputs: full per-model/per-trait GWAS result CSVs from 09.
# Outputs: outputs/GWAS/Combined_Manhattan_plots/<model>/; combined CSV and rMVP JPG plots.
# Paper: Fig. 3a Manhattan and 3b QQ; density/circular plots also retained.
# Analysis section: 10 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_10 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("rMVP"))
  p_thr <- efp_settings$alpha / efp_settings$n_markers
  for (model in efp_settings$models) {
   run_dir <- file.path(project_root, "outputs/GWAS", paste0("results_", model,
                        "_MAF_", efp_settings$maf, "_PCA_total_", efp_settings$pc_total))
   tables <- lapply(efp_settings$traits, function(trait) {
     x <- .efp_results(.efp_result_file(run_dir, model, trait))
     x <- x[, c("SNP", "Chr", "Pos", "P.value")]; names(x)[4] <- trait; x
   })
   merged_gwas <- Reduce(function(a,b) merge(a,b, by = c("SNP", "Chr", "Pos"), all = TRUE), tables)
   chr <- toupper(as.character(merged_gwas$Chr)); is_numeric <- grepl("^[0-9]+$", chr)
   chr[is_numeric] <- wheat_chromosomes[as.integer(chr[is_numeric])]
   if (anyNA(chr) || any(!chr %in% wheat_chromosomes)) stop("Unknown GWAS chromosome code.")
   # rMVP receives numeric chromosomes, with wheat labels in the same order.
   merged_gwas$Chr <- match(chr, wheat_chromosomes)
   merged_gwas <- merged_gwas[order(merged_gwas$Chr, merged_gwas$Pos), ]
   chr_sorted <- wheat_chromosomes[sort(unique(merged_gwas$Chr))]
   folder <- .efp_dir(file.path(project_root, "outputs/GWAS/Combined_Manhattan_plots", model)); setwd(folder)
   write.csv(merged_gwas, paste0("Combined_GWAS_results_", model, ".csv"), row.names = FALSE)
   trait_colors <- c("#0072B2", "#CC0077B3", "#E69F00")
   rMVP::MVP.Report(merged_gwas[, 1:3], plot.type = "d", col = c("darkgreen","yellow","red"), file.type = "jpg", memo = "Density", dpi = 300)
   rMVP::MVP.Report(merged_gwas, plot.type = "q", col = trait_colors, threshold = p_thr,
     signal.pch = 19, signal.cex = 1.5, signal.col = "red", box = FALSE, multracks = TRUE,
     file.type = "jpg", memo = "QQplot", dpi = 300, width = 3, height = 3)
   rMVP::MVP.Report(merged_gwas, plot.type = "m", multracks = TRUE, chr.labels = chr_sorted,
     threshold = c(p_thr, p_thr), threshold.lty = c(1,2), threshold.lwd = c(1,1), threshold.col = c("red","grey"),
     amplify = TRUE, bin.size = 1e6, chr.den.col = c("darkgreen","yellow","red"),
     signal.col = trait_colors, signal.cex = c(1,1), file.type = "jpg", memo = "Manhattan", dpi = 300, width = 14, height = 7)
   rMVP::MVP.Report(merged_gwas, plot.type = "c", r = 0.4, col = c("#082b5e","#1257b8","#99b9e8"),
     chr.labels = chr_sorted, threshold = c(p_thr,p_thr), cir.chr.h = 1.5, amplify = TRUE,
     threshold.lty = c(1,2), threshold.col = c("red","blue"), signal.line = 1,
     signal.col = c("red","green"), chr.den.col = c("darkgreen","yellow","red"), bin.size = 1e6,
     outward = FALSE, file.type = "jpg", memo = "Circular", dpi = 300, width = 5, height = 5)
  }
  invisible(NULL)
}


# =============================================================================
# 11 | Significant markers and sample allele matrix
# Packages: dplyr.
# Inputs: GWAS result CSVs from 09; filtered DB HapMap from 06.
# Outputs: outputs/GWAS/sigMarkers_info.csv, list_of_sigMarkers.csv, transposed_hapmap_with_only_sigMarkers_data.csv.
# Paper: Tables for Fig. 4/6; this selects significant SNPs, not an automated seven-locus clustering.
# Analysis section: 11 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_11 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr"))
  out_dir <- .efp_dir("outputs/GWAS")
  p_thr <- efp_settings$alpha / efp_settings$n_markers
  all_results <- list()
  for (model in efp_settings$models) for (trait in efp_settings$traits) {
   run_dir <- file.path(out_dir, paste0("results_",model,"_MAF_",efp_settings$maf,"_PCA_total_",efp_settings$pc_total))
   x <- .efp_results(.efp_result_file(run_dir, model, trait))
   x$Log10_P.value <- -log10(x$P.value)
   names(x)[names(x) == "H&B.P.Value"] <- "H_B_P.Value"
   x <- x[!is.na(x$P.value) & x$P.value <= p_thr, , drop = FALSE]
   x$Model <- rep(model, nrow(x)); x$Trait <- rep(trait, nrow(x)); names(x)[names(x) == "SNP"] <- "sigMarker"
   all_results[[paste(model,trait)]] <- x
  }
  all_results <- dplyr::bind_rows(all_results)
  write.csv(all_results, file.path(out_dir, "sigMarkers_info.csv"), row.names = FALSE)
  markers <- unique(all_results$sigMarker)
  write.csv(data.frame(sigMarker = markers), file.path(out_dir,"list_of_sigMarkers.csv"), row.names = FALSE)
  if (!length(markers)) { message("No SNPs passed the threshold."); return(invisible(all_results)) }
  h <- .efp_hapmap("inputs/Filtered_EFP_genotypes_AxiomTaNG1_KASP_db_hmp_thresholds_0.1_0.1.txt")
  missing <- setdiff(markers, h[[1]])
  if (length(missing)) stop("Significant markers absent from DB HapMap: ", paste(missing, collapse=", "))
  h <- h[match(markers,h[[1]]), , drop=FALSE]
  samples <- as.data.frame(t(as.matrix(h[, -(1:11), drop=FALSE])), check.names=FALSE, stringsAsFactors=FALSE)
  names(samples) <- markers
  samples <- data.frame(ID=rownames(samples), samples, check.names=FALSE, row.names=NULL)
  write.csv(samples, file.path(out_dir,"transposed_hapmap_with_only_sigMarkers_data.csv"), row.names=FALSE)
  invisible(all_results)
}


# =============================================================================
# 12 | Allele effects and ANOVA
# Packages: dplyr, tidyr, ggplot2, gridExtra, multcompView.
# Inputs: sample-only allele CSV from 11; wide phenotype CSV.
# Outputs: outputs/GWAS/SigMarkers_plots/; three multipage Boxplot*.pdf, merged data, ANOVA means/differences and Tukey CSVs.
# Paper: Fig. 4 panels; select the seven paper markers from multipage PDF. Colours use medians; percent differences use means.
# Analysis section: 12 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_12 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "tidyr", "ggplot2", "gridExtra", "multcompView"))
  ########################################################
  folder <- .efp_dir("outputs/GWAS/SigMarkers_plots")
  SNPdata <- .efp_csv("outputs/GWAS/transposed_hapmap_with_only_sigMarkers_data.csv", "ID")
  Phenodata <- .efp_csv("inputs/EliteFibrePanel_phenotyping.csv", c("ID","RAGT","LIMAGRAIN","Average"))
  .efp_unique_ids(SNPdata$ID, "Allele matrix"); .efp_unique_ids(Phenodata$ID, "Phenotypes")
  variants <- setdiff(names(SNPdata), "ID")
  # Section 11 now retains original punctuation in IDs and marker column names.
  # No '-' to '.' substitutions or deletion of the first ten sample rows.
  for (tr in efp_settings$traits) Phenodata[[tr]] <- .efp_number(Phenodata[[tr]], tr)
  merged_data <- merge(Phenodata[, c("ID", efp_settings$traits)], SNPdata, by="ID", all=FALSE)
  if (!nrow(merged_data)) stop("No matching sample IDs between genotype and phenotype tables.")
  write.csv(data.frame(ID=setdiff(Phenodata$ID,SNPdata$ID)), file.path(folder,"Phenotype_IDs_without_alleles.csv"), row.names=FALSE)
  write.csv(merged_data, file.path(folder,"transposed_hapmap_for_sigMarkers_with_pheyotypic_data.csv"), row.names=FALSE)
  Mydata <- tidyr::pivot_longer(merged_data, cols=all_of(efp_settings$traits), names_to="Trait", values_to="Xylose")
  write.csv(Mydata, file.path(folder,"transposed_long_hapmap_for_sigMarkers_with_pheyotypic_data.csv"), row.names=FALSE)
  set.seed(efp_settings$seed) # reproducible jitter only; statistics unchanged

  # Define output PDF
  fnpdf <- "outputs/GWAS/SigMarkers_plots/Boxplot_SigMarkers_alleles_vs_Xylose.pdf"

  # Define genotype colour mapping
  allele_colors <- c(
    "AA" = "red2",
    "CC" = "blue2",
    "GG" = "green4",
    "TT" = "orange",
    "NN" = "lightgray",
    "AG" = "gray", "GA" = "gray",
    "CT" = "gray", "TC" = "gray",
    "AC" = "gray", "CA" = "gray",
    "GT" = "gray", "TG" = "gray",
    "AT" = "gray", "TA" = "gray",
    "CG" = "gray", "GC" = "gray"
  )

  # Collect all plots in a list
  plot_list <- list()

  for (variant in variants) {
    p <- ggplot(Mydata, aes(x = .data[[variant]], y = Xylose, fill = .data[[variant]])) +
      geom_boxplot() +
      geom_jitter(width = 0.1, height = 0, alpha = 0.3, shape = 21) +
      facet_wrap(~Trait, nrow = 1) +
      ggtitle(variant) +
      xlab("Allele") +
      ylab("Xylose equivalents (mg/g)") +
      scale_fill_manual(values = allele_colors, na.value = "gray90") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "none"
      )
    plot_list[[variant]] <- p
  }
  # Arrange 4 plots per page and save to PDF
  .efp_multipage(plot_list, fnpdf, width=4, height=4)

  # --------------------------------------------------
  # Boxplots homozygous alleles only
  # --------------------------------------------------
  # Output PDF
  fnpdf <- "outputs/GWAS/SigMarkers_plots/Boxplot_SigMarkers_alleles_vs_Xylose_homozygous_alleles_only.pdf"

  # Collect plots
  plot_list <- list()

  # Valid homozygous genotypes (keep your set)
  homozygous_genotypes <- c("AA", "CC", "GG", "TT")

  for (variant in variants) {
  
    # Filter to keep only homozygous alleles
    filtered_data <- Mydata %>%
      dplyr::filter(.data[[variant]] %in% homozygous_genotypes) %>%
      dplyr::mutate(Allele = .data[[variant]])
  
    if (nrow(filtered_data) == 0) next  # Skip if no valid data
  
    # Compute medians per Trait x Allele
    med_tbl <- filtered_data %>%
      dplyr::group_by(Trait, Allele) %>%
      dplyr::summarise(med = stats::median(Xylose, na.rm = TRUE), .groups = "drop")
  
    # For each Trait: mark highest/lowest allele by median
    med_tbl <- med_tbl %>%
      dplyr::group_by(Trait) %>%
      dplyr::mutate(
        rank_med = dplyr::dense_rank(med),           # 1 = lowest
        rank_med_rev = dplyr::dense_rank(-med),      # 1 = highest
        Colour = dplyr::case_when(
          rank_med_rev == 1 ~ "#68228B",  # highest median
          rank_med == 1     ~ "#d18df2",  # lowest median
          TRUE              ~ "grey80"    # middle (if >2 alleles exist)
        )
      ) %>%
      dplyr::ungroup()
  
    # Join colours back to the plotting data (Trait-specific colouring)
    filtered_data2 <- filtered_data %>%
      dplyr::left_join(med_tbl %>% dplyr::select(Trait, Allele, Colour),
                       by = c("Trait", "Allele"))
  
    # Use identity scale because Colour column already contains hex codes
    p <- ggplot(filtered_data2, aes(x = Allele, y = Xylose, fill = Colour)) +
      geom_boxplot() +
      geom_jitter(width = 0.1, height = 0, alpha = 0.3, shape = 21) +
      facet_wrap(~Trait, nrow = 1) +
      ggtitle(variant) +
      xlab("Allele") +
      ylab("Xylose equivalents (mg/g)") +
      scale_fill_identity() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "none"
      )
  
    plot_list[[variant]] <- p
  }

  # Plot 1 per page (as in your code)
  .efp_multipage(plot_list, fnpdf, width=4, height=4)

  ####################################################
  # Boxplots homozygous alleles with stats (Fig. 4)
  ####################################################

  # --------------------------------------------------
  # Output files
  # --------------------------------------------------
  fnpdf <- "outputs/GWAS/SigMarkers_plots/Boxplot_SigMarkers_alleles_vs_Xylose_homozygous_alleles_only_with_stats.pdf"
  fncsv <- "outputs/GWAS/SigMarkers_plots/Allele_stats_ANOVA_meanDiff_homozygous_only.csv"

  dir.create(dirname(fnpdf), recursive = TRUE, showWarnings = FALSE)

  plot_list <- list()
  all_stats  <- list()
  all_tukey <- list()

  homozygous_genotypes <- c("AA", "CC", "GG", "TT")

  # --------------------------------------------------
  # Loop over variants
  # --------------------------------------------------
  for (variant in variants) {
  
    filtered_data <- Mydata %>%
      filter(.data[[variant]] %in% homozygous_genotypes) %>%
      mutate(Allele = .data[[variant]]) %>%
      filter(!is.na(Xylose), !is.na(Trait))
  
    if (nrow(filtered_data) == 0) next
  
    # ------------------------------------------------
    # Median-based colouring
    # ------------------------------------------------
    med_tbl <- filtered_data %>%
      group_by(Trait, Allele) %>%
      summarise(med = median(Xylose, na.rm = TRUE), .groups = "drop") %>%
      group_by(Trait) %>%
      mutate(
        Colour = case_when(
          med == max(med) ~ "#68228B",
          med == min(med) ~ "#d18df2",
          TRUE ~ "grey80"
        )
      ) %>%
      ungroup()
  
    plot_data <- filtered_data %>%
      left_join(med_tbl %>% select(Trait, Allele, Colour),
                by = c("Trait", "Allele"))
  
    # ------------------------------------------------
    # Fixed annotation bands (avoid overlap)
    # ------------------------------------------------
    y_limits_tbl <- plot_data %>%
      group_by(Trait) %>%
      summarise(
        y_max = max(Xylose, na.rm = TRUE),
        y_min = min(Xylose, na.rm = TRUE),
        y_rng = y_max - y_min,
        y_p   = y_max + 0.18 * y_rng,
        y_lab = y_max + 0.08 * y_rng,
        .groups = "drop"
      )
  
    # ------------------------------------------------
    # Mean difference + ANOVA + letters
    # ------------------------------------------------
    mean_tbl <- plot_data %>%
      group_by(Trait, Allele) %>%
      summarise(
        n = n(),
        mean_xylose = mean(Xylose, na.rm = TRUE),
        .groups = "drop"
      )
  
    trait_summary <- mean_tbl %>%
      group_by(Trait) %>%
      summarise(
        allele_low  = Allele[which.min(mean_xylose)][1],
        allele_high = Allele[which.max(mean_xylose)][1],
        mean_low  = min(mean_xylose),
        mean_high = max(mean_xylose),
        mean_diff_percent = (mean_high - mean_low) / mean_low * 100,
        .groups = "drop"
      )
  
    letters_all <- NULL
    anova_rows  <- list()
  
    for (tr in unique(plot_data$Trait)) {
    
      sub_df <- plot_data %>% filter(Trait == tr)
      sub_df$Allele <- droplevels(factor(sub_df$Allele))
    
      pval <- NA_real_
    
      if (nlevels(sub_df$Allele) >= 2 && nrow(sub_df) > nlevels(sub_df$Allele) && sd(sub_df$Xylose) > 0) {
        aov_fit <- aov(Xylose ~ Allele, data = sub_df)
        pval <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
      
        tuk <- TukeyHSD(aov_fit, "Allele")
        tuk_p <- tuk$Allele[, "p adj"]
        names(tuk_p) <- rownames(tuk$Allele)
        all_tukey[[paste(variant,tr)]] <- data.frame(Variant=variant, Trait=tr, Comparison=names(tuk_p), Adjusted_P=as.numeric(tuk_p))
      
        cld <- multcompView::multcompLetters(tuk_p, threshold = 0.05)
      
        letters_df <- data.frame(
          Trait  = tr,
          Allele = names(cld$Letters),
          label  = cld$Letters,
          stringsAsFactors = FALSE
        ) %>%
          left_join(y_limits_tbl %>% select(Trait, y_lab), by = "Trait") %>%
          rename(y = y_lab)
      
        letters_all <- bind_rows(letters_all, letters_df)
      }
    
      anova_rows[[tr]] <- data.frame(Trait = tr, p_value = pval)
    }
  
    anova_df <- bind_rows(anova_rows)
  
    stats_out <- trait_summary %>%
      left_join(anova_df, by = "Trait") %>%
      mutate(Variant = variant)
  
    all_stats[[variant]] <- stats_out
  
    # ------------------------------------------------
    # p-value text (fixed band)
    # ------------------------------------------------
    pval_df <- anova_df %>%
      mutate(
        p_text = case_when(
          is.na(p_value) ~ "p = NA",
          p_value < 0.001 ~ "p < 0.001",
          TRUE ~ paste0("p = ", signif(p_value, 3))
        )
      ) %>%
      left_join(y_limits_tbl %>% select(Trait, y_p), by = "Trait")
  
    # ------------------------------------------------
    # Plot
    # ------------------------------------------------
    p <- ggplot(plot_data, aes(x = Allele, y = Xylose, fill = Colour)) +
      geom_boxplot() +
      geom_jitter(width = 0.1, alpha = 0.3, shape = 21) +
      facet_wrap(~Trait, nrow = 1) +
      scale_fill_identity() +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      labs(
        title = variant,
        x = "Allele",
        y = "Xylose equivalents (mg/g)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "none"
      )
  
    if (!is.null(letters_all)) {
      p <- p +
        geom_text(
          data = letters_all,
          aes(x = Allele, y = y, label = label),
          size = 5,
          vjust = 0,
          inherit.aes = FALSE
        )
    }
  
    p <- p +
      geom_text(
        data = pval_df,
        aes(x = -Inf, y = y_p, label = p_text),
        hjust = -0.05,
        vjust = 0,
        size = 4,
        inherit.aes = FALSE
      )
  
    plot_list[[variant]] <- p
  }

  # --------------------------------------------------
  # Save outputs
  # --------------------------------------------------
  stats_df <- bind_rows(all_stats)
  write.csv(stats_df, fncsv, row.names = FALSE)
  write.csv(bind_rows(all_tukey), file.path(folder,"Allele_Tukey_comparisons.csv"), row.names=FALSE)

  .efp_multipage(plot_list, fnpdf, width=4, height=4)
  invisible(stats_df)
}


# =============================================================================
# 13 | High-fibre allele frequencies
# Packages: dplyr, tidyr, ggplot2.
# Inputs: curated 1_SigSNPs_to_calculate_allele_frequency.csv and 4_SigMarkers_allele_frequency_for_plotting.csv in outputs/GWAS/AlleleFrequency/.
# Outputs: same folder: 2_Allele_counts_per_group_per_SNP.csv, allele_frequencies_in_EFP_for_7_marker.pdf, allele_frequencies_by_7_marker.pdf.
# Paper: Fig. 5a and S4; HFA/LFA assignments and the percentage table are supplied inputs.
# Analysis section: 13 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_13 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "tidyr", "ggplot2"))
  # Table 4 is NOT generated from table 2 by the supplied code.
  # HFA/LFA frequencies in table 4 must exclude heterozygous and missing calls.
  df <- .efp_csv("outputs/GWAS/AlleleFrequency/1_SigSNPs_to_calculate_allele_frequency.csv", c("ID", "group"))
  # Columns that contain SNP calls
  snp_cols <- setdiff(colnames(df), c("ID", "group"))

  # Convert to long format
  long_df <- df %>%
    pivot_longer(
      cols = all_of(snp_cols),
      names_to = "SNP_name",
      values_to = "Call"
    )

  # Count HFA / LFA / Het / NN per group and per SNP
  counts_df <- long_df %>%
    group_by(group, SNP_name, Call) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = Call,
      values_from = Count,
      values_fill = 0
    ) %>%
    arrange(SNP_name, group)

  # View result
  print(counts_df)

  # Save to CSV
  write.csv(
    counts_df,
    "outputs/GWAS/AlleleFrequency/2_Allele_counts_per_group_per_SNP.csv",
    row.names = FALSE
  )

  #################################################################################
  # High-Fibre Allele frequencies in EFP only (Fig. 5a)
  #---------------------------------------------------------------
  # Load libraries

  # Set the working directory as the parent directory of the R script's folder
  setwd(project_root)

  # ==== Step 1: Read the dataset from CSV ====
  df <- .efp_csv("outputs/GWAS/AlleleFrequency/4_SigMarkers_allele_frequency_for_plotting.csv", c("Marker","GeneticResource","Allele","AlleleFrequency"))
  df$AlleleFrequency <- .efp_number(df$AlleleFrequency, "AlleleFrequency")
  if (any(df$AlleleFrequency < 0 | df$AlleleFrequency > 100, na.rm=TRUE)) stop("AlleleFrequency must be percent (0-100).")

  # ==== Step 2: Filter + clean (keep EFP only) ====
  df_clean <- df %>%
    filter(
      !is.na(AlleleFrequency),
      GeneticResource == "EFP"
    ) %>%
    mutate(
      Marker = factor(Marker, levels = rev(unique(Marker))),  # <-- alphabetical order
      Allele = factor(Allele, levels = c("HFA", "LFA"))
    )

  # ==== Step 3: Create the stacked HORIZONTAL bar plot ====
  p <- ggplot(df_clean, aes(x = AlleleFrequency, y = Marker, fill = Allele)) +
    geom_col(width = 0.9, position = position_stack(reverse = TRUE)) +
    #facet_wrap(~ Marker, scales = "fixed", nrow = 2) +
    scale_fill_manual(
      values = c(
        "HFA" = "#68228B",
        "LFA" = "#d18df2"
      ),
      drop = FALSE
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.text.x = element_text(),
      panel.grid.major.y = element_blank()
    ) +
    labs(
      title = "Allele Frequencies per Marker",
      x = "Allele Frequency (%)",
      y = "Significant SNP",
      fill = "Allele"
    )

  # ==== Step 4: Save to PDF ====
  .efp_ggsave("outputs/GWAS/AlleleFrequency/allele_frequencies_in_EFP_for_7_marker.pdf", plot = p, width = 6, height = 3)

  #####################################################################
  # High-Fibre Allele frequencies in EFP, Variety, Landrace (Fig. S4)
  #====================================================================

  # Load libraries

  # Set the working directory as the parent directory of the R script's folder
  setwd(project_root)

  # ==== Step 1: Read the dataset from CSV ====
  df <- .efp_csv("outputs/GWAS/AlleleFrequency/4_SigMarkers_allele_frequency_for_plotting.csv", c("Marker","GeneticResource","Allele","AlleleFrequency"))
  df$AlleleFrequency <- .efp_number(df$AlleleFrequency, "AlleleFrequency")
  if (any(df$AlleleFrequency < 0 | df$AlleleFrequency > 100, na.rm=TRUE)) stop("AlleleFrequency must be percent (0-100).")
  # ==== Step 2: Clean + set order (EFP, Variety, Landrace) ====
  df_clean <- df %>%
    filter(!is.na(AlleleFrequency)) %>%
    mutate(
      GeneticResource = factor(GeneticResource, levels = c("EFP", "Variety", "Landrace")),
      Allele = factor(Allele, levels = c("HFA", "LFA"))  # optional: consistent fill order
    )

  # ==== Step 3: Create the stacked HORIZONTAL bar plot ====
  p <- ggplot(df_clean, aes(x = AlleleFrequency, y = GeneticResource, fill = Allele)) +
    geom_col(width = 0.9) +
    facet_wrap(~ Marker, scales = "fixed", nrow = 2) +
    scale_fill_manual(
      values = c(
        "HFA" = "#68228B",
        "LFA" = "#d18df2"
      ),
      drop = FALSE
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.text.x = element_text(),
      panel.grid.major.y = element_blank()
    ) +
    labs(
      title = "Allele Frequencies per Marker",
      x = "Allele Frequency (%)",
      y = "Genetic Resource",
      fill = "Allele"
    )

  # ==== Step 4: Save to PDF ====
  .efp_ggsave("outputs/GWAS/AlleleFrequency/allele_frequencies_by_7_marker.pdf", plot = p, width = 6, height = 3)
}


# =============================================================================
# 14 | Favourable allele counts and phenotypes
# Packages: dplyr, tidyr, readr, tibble, ggplot2, multcompView.
# Inputs: curated 5_for_HFA_counting_per_EFP_line.csv and 7_HFA_counts_per_EFP_line_with_average_Xylose.csv in AlleleFrequency/.
# Outputs: same folder: 6_HFA_counts_per_EFP_line.csv, Boxplot_Average_Xylose_by_HFA_count_7_markers.pdf, ANOVA/Tukey CSVs.
# Paper: Fig. 5b; ANOVA across count classes, Tukey adjusted alpha=0.05 (not 0.001).
# Analysis section: 14 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_14 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "tidyr", "readr", "tibble", "ggplot2", "multcompView"))
  set.seed(efp_settings$seed)
  # Table 7 is a separate supplied phenotype join; table 6 alone does not create it.
  # This is an observational comparison, not a breeding-response prediction.
  # ==== Step 1: Read the dataset from CSV ====
  df <- .efp_csv("outputs/GWAS/AlleleFrequency/5_for_HFA_counting_per_EFP_line.csv", "ID")

  # Identify SNP columns (everything except ID)
  snp_cols <- setdiff(names(df), "ID")

  # Count number of HFA per row
  df$N_missing_calls <- rowSums(sapply(df[,snp_cols,drop=FALSE], function(x) is.na(x) | x %in% c("NN", "NA", "")))
  # Audit missing calls; the historical HFA count is a lower bound for incomplete rows.
  df_hfa_counts <- df %>%
    mutate(
      HFA_count = rowSums(across(
        all_of(snp_cols),
        ~ .x == "HFA"
      ), na.rm = TRUE)
    )

  # View result

  # Optionally write to file
  write.csv(df_hfa_counts,
            "outputs/GWAS/AlleleFrequency/6_HFA_counts_per_EFP_line.csv",
            row.names = FALSE)

  #=========================================================
  # plot Average Xylose values against the number of HFAs
  #=========================================================

  # Set the working directory as the parent directory of the R script's folder
  setwd(project_root)

  # -----------------------------
  # Input / output
  # -----------------------------
  in_csv <- "outputs/GWAS/AlleleFrequency/7_HFA_counts_per_EFP_line_with_average_Xylose.csv"  # <-- change this
  out_pdf <- "outputs/GWAS/AlleleFrequency/Boxplot_Average_Xylose_by_HFA_count_7_markers.pdf"
  out_csv <- "outputs/GWAS/AlleleFrequency/ANOVA_Tukey_letters_Average_Xylose_by_HFA_count.csv"

  # -----------------------------
  # Read + robust clean
  # -----------------------------
  df_raw <- .efp_csv(in_csv, c("HFA_count","Average_Xylose"))

  # Make sure the required columns exist (exact names)
  req_cols <- c("HFA_count", "Average_Xylose")
  missing_cols <- setdiff(req_cols, colnames(df_raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "))
  }

  df <- df_raw %>%
    mutate(
      # robust numeric parsing:
      # works for "5", "5 HFAs", "HFA=5", " 5 ", etc.
      HFA_count_num = .efp_number(HFA_count, "HFA_count"),
      Average_Xylose = .efp_number(Average_Xylose, "Average_Xylose")
    ) %>%
    filter(!is.na(HFA_count_num), !is.na(Average_Xylose)) %>%
    filter(HFA_count_num %in% 0:6) %>%
    mutate(
      HFA_count_f = factor(
        HFA_count_num,
        levels = 6:0,
        labels = paste0(6:0, " HFAs")
      )
    )

  if (any(!is.na(df_raw$HFA_count) & !df_raw$HFA_count %in% 0:6))
    stop("The paper plotting code expects HFA counts 0-6; revise explicitly for another dataset.")
  # If empty, fail early with helpful message
  if (nrow(df) == 0) {
    stop(
      "After cleaning/filtering, df has 0 rows.\n",
      "Check that HFA_count contains values 0–6 and Average_Xylose is numeric.\n",
      "Tip: print(unique(df_raw$HFA_count)) and head(df_raw$Average_Xylose)."
    )
  }

  # -----------------------------
  # Gradient palette (6 -> 0)
  # -----------------------------
  pal <- colorRampPalette(c("#4d0273", "#eedcf5"))(7)
  names(pal) <- paste0(6:0, " HFAs")

  # -----------------------------
  # ANOVA + Tukey + letters
  # -----------------------------
  anova_p <- NA_real_
  letters_df <- NULL

  if (nlevels(droplevels(df$HFA_count_f)) >= 2 && length(unique(df$Average_Xylose)) > 1) {
    aov_fit <- aov(Average_Xylose ~ HFA_count_f, data = df)
    anova_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  
    tuk <- TukeyHSD(aov_fit, "HFA_count_f")
    tuk_p <- tuk$HFA_count_f[, "p adj"]
    names(tuk_p) <- rownames(tuk$HFA_count_f)
    write.csv(data.frame(Comparison=names(tuk_p), Adjusted_P=as.numeric(tuk_p)), "outputs/GWAS/AlleleFrequency/HFA_Tukey_comparisons.csv", row.names=FALSE)
  
    cld <- multcompView::multcompLetters(tuk_p, threshold = 0.05)
  
    # y position for letters above each box
    y_max <- tapply(df$Average_Xylose, df$HFA_count_f, max, na.rm = TRUE)
    rng <- range(df$Average_Xylose, na.rm = TRUE)
    offset <- 0.06 * (rng[2] - rng[1])
  
    letters_df <- tibble(
      HFA_count_f = factor(names(cld$Letters), levels = levels(df$HFA_count_f)),
      letter = cld$Letters,
      y = unname(y_max[names(cld$Letters)]) + offset
    )
  }

  # -----------------------------
  # Save stats table
  # -----------------------------
  stats_tbl <- df %>%
    group_by(HFA_count_f) %>%
    summarise(
      n = n(),
      mean_Average_Xylose = mean(Average_Xylose, na.rm = TRUE),
      sd_Average_Xylose   = sd(Average_Xylose, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      anova_p_value = anova_p,
      anova_p_text = case_when(
        is.na(anova_p) ~ "p = NA",
        anova_p < 0.001 ~ "p < 0.001",
        TRUE ~ paste0("p = ", signif(anova_p, 3))
      )
    ) %>%
    left_join(
      if (!is.null(letters_df)) letters_df %>% select(HFA_count_f, letter)
      else tibble(HFA_count_f = levels(df$HFA_count_f), letter = NA_character_),
      by = "HFA_count_f"
    )

  write.csv(stats_tbl, out_csv, row.names = FALSE)

  # -----------------------------
  # Plot (p-value top-left; letters above boxes)
  # -----------------------------
  p_text <- case_when(
    is.na(anova_p) ~ "ANOVA: p = NA",
    anova_p < 0.001 ~ "ANOVA: p < 0.001",
    TRUE ~ paste0("ANOVA: p = ", signif(anova_p, 3))
  )

  p <- ggplot(df, aes(x = HFA_count_f, y = Average_Xylose, fill = HFA_count_f)) +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.3) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1, shape = 21, stroke = 0) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(
      x = "Number of high-fibre alleles per line",
      y = "Average Xylose equivalents (mg/g)",
      title = "Average xylose by number of high-fibre alleles"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.line = element_line(color = "black")
    ) +
    annotate("text", x = -Inf, y = Inf, label = p_text, hjust = -0.05, vjust = 1.2, size = 4)

  if (!is.null(letters_df)) {
    p <- p +
      geom_text(
        data = letters_df,
        aes(x = HFA_count_f, y = y, label = letter),
        inherit.aes = FALSE,
        size = 5,
        vjust = 0
      )
  }

  .efp_ggsave(out_pdf, p, width = 4, height = 4)
  print(p)

  cat("Rows plotted:", nrow(df), "\n")
  cat("Saved:\n", out_pdf, "\n", out_csv, "\n")
}


# =============================================================================
# 15 | Chromosome LD plots matching the supplied supplementary layout
# Packages: snpStats, ggplot2, scales.
# Inputs: numerical genotype/map from 07a.
# Outputs: outputs/LD_full_chromosome/LD_pairs_chr*.csv, LD_heatmap_chr*.pdf, LD_triangular_R2_chr*.{pdf,jpg}, LD_decay_chr*.pdf and bin CSVs.
# Paper: S6 physical-position heatmaps and S7 original adaptive bins; see 16 for SNP-spaced alternatives.
# Analysis section: 15 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_15 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("snpStats", "ggplot2", "scales"))
  ld_data <- read_ld_data("outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.Numerical.txt",
                          "outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.map.txt")
  # Original bin width starts at max(5 x median marker spacing, 50 kb).
  # The target 30 is the median occupied-bin count; sparse bins remain plotted.
  # Reuses the latest pairwise validity checks; undefined LD remains NA.
  # ---- Output folder
  out_dir <- "outputs/LD_full_chromosome"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---------- Helper: robust LD-decay summary using numeric bin indices ----------
  decay_summary_numeric_bins <- function(pairs_df,
                                         target_pairs_per_bin = 30L,
                                         initial_bin_bp = 5e4,
                                         max_bin_bp = NULL) {
    if (nrow(pairs_df) == 0) return(NULL)
    df <- subset(pairs_df, is.finite(Distance_bp) & is.finite(R2) & Distance_bp >= 0)
    if (nrow(df) == 0) return(NULL)
  
    if (is.null(max_bin_bp)) {
      span <- max(df$Distance_bp, na.rm = TRUE)
      max_bin_bp <- max(1e5, 0.02 * span)
    }
  
    bin_width_bp <- max(1e3, initial_bin_bp)
  
    compute_stats <- function(df, bw) {
      if (nrow(df) == 0) return(NULL)
      df$bin_index <- floor(df$Distance_bp / bw)
      means  <- aggregate(R2 ~ bin_index, df, mean)
      counts <- aggregate(R2 ~ bin_index, df, length)
      colnames(counts)[2] <- "Pair_Count"
      out <- merge(means, counts, by = "bin_index")
      out$Distance_mid_Mb <- ((out$bin_index + 0.5) * bw) / 1e6
      out
    }
  
    repeat {
      stat <- compute_stats(df, bin_width_bp)
      if (is.null(stat) || nrow(stat) == 0) break
      med_pairs <- stats::median(stat$Pair_Count, na.rm = TRUE)
      if (is.finite(med_pairs) && med_pairs >= target_pairs_per_bin) break
      new_bw <- ceiling(bin_width_bp * 1.5)
      if (new_bw >= max_bin_bp) {
        bin_width_bp <- max_bin_bp
        stat <- compute_stats(df, bin_width_bp)
        break
      }
      bin_width_bp <- new_bw
    }
  
    list(summary = stat, bin_width_bp = bin_width_bp)
  }

  # ---------- Helper: high-contrast LD heatmap with R^2 legend (margin-safe) ----------
  # ---------- Mb-axis LD heatmap (equal x/y scale; irregular grid safe) ----------
  plot_ld_heatmap <- function(ld_mat, snps_chr, chr, out_pdf_heat, chr_info) {
    pal <- colorRampPalette(c("#FFFFFF", "#FFE082", "#FB8C00", "#E53935", "#7F0000"))
  
    # order by SNP order; keep only one triangle
    M <- suppressWarnings(as.matrix(ld_mat))
    M <- M[snps_chr, snps_chr, drop = FALSE]
    M[lower.tri(M, diag = TRUE)] <- NA_real_   # show upper triangle only
    M[!is.finite(M)] <- NA_real_
    n <- nrow(M); if (n < 2) { message("Chr ", chr, ": too few SNPs for heatmap."); return(invisible()) }
  
    # Mb positions (ensure strictly increasing even if there are ties)
    pos_mb <- as.numeric(chr_info$Position[match(snps_chr, chr_info$SNP)]) / 1e6
    if (anyDuplicated(pos_mb)) {
      pos_mb <- pos_mb + ave(pos_mb, pos_mb, FUN = function(v) (seq_along(v) - 1) * 1e-9)
    }
    rng <- range(pos_mb, finite = TRUE)
  
    # nice Mb ticks
    at_mb <- pretty(rng, n = 8); at_mb <- at_mb[at_mb >= rng[1] & at_mb <= rng[2]]
  
    pdf(out_pdf_heat, width = 7.5, height = 6, bg = "white")
    old_par <- par(no.readonly = TRUE); on.exit({ try(par(old_par), silent = TRUE); dev.off() }, add = TRUE)
  
    # main panel — plot on Mb grid; equal scaling; smaller Mb at top
    par(fig = c(0.00, 0.86, 0.00, 1.00), mar = c(4.2, 4.6, 3.0, 0.2), xaxs = "i", yaxs = "i")
    image(x = pos_mb, y = pos_mb, z = t(M),             # transpose matches expected orientation
          col = pal(256), zlim = c(0, 1),
          axes = FALSE, xlab = "Physical position (Mb)", ylab = "Physical position (Mb)",
          asp = 1, ylim = rev(rng))                    # equal units; flip y so triangle points down
    axis(1, at = at_mb, labels = sprintf("%.2f", at_mb), las = 1, cex.axis = 0.9)
    axis(2, at = at_mb, labels = sprintf("%.2f", at_mb), las = 1, cex.axis = 0.9)
    box()
    title(main = bquote("LD ("*R^2*") heatmap — chr "*.(chr)))
  
    # legend
    par(fig = c(0.88, 0.98, 0.15, 0.85), mar = c(0.2, 0.2, 0.2, 2.2), new = TRUE, xaxs = "i", yaxs = "i")
    legend_x <- c(0, 1); legend_y <- seq(0, 1, length.out = 257)
    legend_z <- matrix(seq(0, 1, length.out = 256), nrow = 1)
    image(x = legend_x, y = legend_y, z = legend_z, col = pal(256), zlim = c(0, 1), axes = FALSE)
    axis(4, at = seq(0, 1, by = 0.2), labels = sprintf("%.2f", seq(0, 1, by = 0.2)), las = 1, cex.axis = 0.8)
    mtext(expression(R^2), side = 4, line = 1.2, cex = 0.9)
    box()
  
    message("Saved: ", out_pdf_heat)
  }

  # ---------- Helper: triangular R² LD plot (upside-down, panel-aligned Mb axis) ----------
  # Saves both PDF (out_pdf) and a JPG with the same name (.jpg extension).
  # Draws the Mb axis manually so it renders identically in PDF/JPG.
  plot_ld_triangle_R2 <- function(ld_mat, snps_chr, chr, chr_info, out_pdf) {
    pal <- colorRampPalette(c("#FFFFFF", "#FFE082", "#FB8C00", "#E53935", "#7F0000"))
  
    # Order matrix by chromosome SNP order
    M <- suppressWarnings(as.matrix(ld_mat))
    M <- M[snps_chr, snps_chr, drop = FALSE]
    n <- nrow(M)
    if (n < 2) { message("Chr ", chr, ": too few SNPs for triangular plot."); return(invisible()) }
  
    # Upper triangle only
    idx <- which(upper.tri(M), arr.ind = TRUE)
    df <- data.frame(i = idx[,1], j = idx[,2], r2 = as.numeric(M[idx]), stringsAsFactors = FALSE)
    df <- df[is.finite(df$r2), ]
    if (!nrow(df)) { message("Chr ", chr, ": no finite R² values to plot."); return(invisible()) }
  
    # 45° rotation coords, then flip vertically so the triangle points DOWN
    df$x <- (df$i + df$j) / 2
    df$y <- (df$j - df$i)
    y_max <- (n - 1)
    df$y_flip <- y_max - df$y + 1
  
    # Mb positions for the custom top axis (use pretty breaks, then snap to nearest SNP)
    pos_vec <- chr_info$Position[match(snps_chr, chr_info$SNP)]
    pos_mb  <- pos_vec / 1e6
    pretty_mb <- pretty(range(pos_mb, finite = TRUE), n = 8)
    pretty_mb <- pretty_mb[pretty_mb >= min(pos_mb) & pretty_mb <= max(pos_mb)]
    # nearest SNP index for each pretty Mb value
    tick_idx <- vapply(pretty_mb, function(v) which.min(abs(pos_mb - v)), integer(1))
    # ensure strictly increasing & unique
    o <- order(tick_idx); tick_idx <- unique(tick_idx[o]); tick_mb <- pos_mb[tick_idx]
    tick_lab <- sprintf("%.2f", tick_mb)
  
    # Plot: bigger tiles to make haploblocks clearer
    p <- ggplot(df, aes(x = x, y = y_flip, fill = r2)) +
      geom_tile(width = 2, height = 2) +
      scale_fill_gradientn(colours = pal(256), limits = c(0, 1),
                           name = expression(R^2), na.value = NA) +
      coord_fixed() +
      scale_x_continuous(limits = c(1, n), breaks = tick_idx, labels = NULL, expand = c(0,0)) +
      scale_y_continuous(limits = c(1, y_max), expand = c(0,0)) +
      labs(title = paste0("LD (R^2) triangular plot — chr ", chr),
           x = "SNP order (by position)", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(36, 16, 28, 16),   # a touch more top room
            plot.title  = element_text(margin = margin(b = 10)))
  
    draw_with_axis <- function(device_fun, file, width_in, height_in, res = NULL, bg = "white") {
      if (is.null(res)) device_fun(file, width = width_in, height = height_in, bg = bg) else
        device_fun(file, width = width_in, height = height_in, units = "in", res = res, bg = bg)
    
      g <- ggplotGrob(p)
      grid::grid.newpage(); grid::grid.draw(g)
    
      # Locate panel and open a small viewport exactly above its top edge
      pan   <- g$layout[g$layout$name == "panel", ][1, ]
      left  <- grid::convertWidth(sum(g$widths[1:(pan$l - 1)]), "npc", valueOnly = TRUE)
      width <- grid::convertWidth(sum(g$widths[pan$l:pan$r]),   "npc", valueOnly = TRUE)
      top   <- grid::convertHeight(sum(g$heights[1:(pan$t - 1)]), "npc", valueOnly = TRUE)
    
      grid::pushViewport(
        grid::viewport(x = left, y = 1 - top, width = width,
                       height = grid::unit(2.4, "lines"), just = c("left","top"),
                       gp = grid::gpar(col = "black"))
      )
    
      # Baseline
      grid::grid.lines(x = grid::unit(c(0, 1), "npc"),
                       y = grid::unit(0.4, "lines"))
    
      # Map SNP-index ticks to [0,1] along the top edge
      at <- (tick_idx - 1) / (n - 1)
    
      # Tick marks + labels (draw one-by-one for raster robustness)
      for (i in seq_along(at)) {
        grid::grid.segments(x0 = grid::unit(at[i], "npc"), x1 = grid::unit(at[i], "npc"),
                            y0 = grid::unit(0.4, "lines"), y1 = grid::unit(1.0, "lines"))
        grid::grid.text(label = tick_lab[i],
                        x = grid::unit(at[i], "npc"),
                        y = grid::unit(1.6, "lines"),
                        just = c("center", "bottom"),
                        gp = grid::gpar(cex = 0.9))
      }
    
      # Axis title (centered, above labels)
      grid::grid.text("Physical position (Mb)",
                      x = grid::unit(0.5, "npc"),
                      y = grid::unit(2.2, "lines"),
                      just = c("center", "bottom"),
                      gp = grid::gpar(cex = 0.95))
    
      grid::popViewport()
      dev.off()
    }
  
    # Save PDF
    draw_with_axis(grDevices::pdf, out_pdf, width_in = 10, height_in = 6, bg = "white")
    message("Saved: ", out_pdf)
  
    # Save JPG (same name with .jpg)
    out_jpg <- sub("\\.pdf$", ".jpg", out_pdf, ignore.case = TRUE)
    draw_with_axis(grDevices::jpeg, out_jpg, width_in = 10, height_in = 6, res = 300, bg = "white")
    message("Saved: ", out_jpg)
  }

  # ---- Parameters
  MAX_SNPS_FULL <- 1500L
  TARGET_PAIRS_PER_BIN <- 30L

  # ---- Main loop per chromosome
  chroms <- unique(ld_data$map$Chromosome)
  for (chr in chroms) {
    # Use the common validated sample-by-SNP calculation (no transposition).
    ld_object <- calculate_ld(ld_data, chr, max_snps = MAX_SNPS_FULL)
    chr_info <- ld_object$map; snps_chr <- chr_info$SNP; n_snps <- length(snps_chr)
    if (n_snps < 2L) next
    ld_mat <- ld_object$r2
    # ---- Save pairwise R² (upper triangle) to CSV
    ut <- which(upper.tri(ld_mat), arr.ind = TRUE)
    if (nrow(ut) == 0) { message("Chr ", chr, ": too few SNPs to extract LD; skipping."); next }
  
    snp1 <- rownames(ld_mat)[ut[,1]]
    snp2 <- colnames(ld_mat)[ut[,2]]
    r2   <- ld_mat[ut]
  
    pos_map <- setNames(chr_info$Position, chr_info$SNP)
    pairs_df <- data.frame(
      Chromosome = chr,
      SNP1 = snp1,
      SNP2 = snp2,
      Position1 = as.numeric(pos_map[snp1]),
      Position2 = as.numeric(pos_map[snp2]),
      Distance_bp = abs(as.numeric(pos_map[snp2]) - as.numeric(pos_map[snp1])),
      R2 = as.numeric(r2),
      stringsAsFactors = FALSE
    )
    pairs_df <- pairs_df[is.finite(pairs_df$R2) & !is.na(pairs_df$Position1) & !is.na(pairs_df$Position2), ]
  
    out_csv <- file.path(out_dir, sprintf("LD_pairs_chr%s.csv", chr))
    write.csv(pairs_df, out_csv, row.names = FALSE)
    message("Saved: ", out_csv, " (", nrow(pairs_df), " pairs)")
  
    # ---- Produce heatmap PDF (this was missing before)
    out_pdf_heat <- file.path(out_dir, sprintf("LD_heatmap_chr%s.pdf", chr))
    plot_ld_heatmap(ld_mat, snps_chr, chr, out_pdf_heat, chr_info)
  
  
    # Triangular R² plot (Haploview-like)
    out_pdf_tri <- file.path(out_dir, sprintf("LD_triangular_R2_chr%s.pdf", chr))
    plot_ld_triangle_R2(ld_mat, snps_chr, chr, chr_info, out_pdf_tri)
  
  
    # ---- (Optional) LD decay — keep commented if not needed
    med_spacing <- if (n_snps >= 3) median(diff(sort(chr_info$Position)), na.rm = TRUE) else NA
    init_bin <- if (is.finite(med_spacing)) max(5 * med_spacing, 5e4) else 5e4
    decay <- decay_summary_numeric_bins(pairs_df, TARGET_PAIRS_PER_BIN, init_bin)
    if (!is.null(decay$summary) && nrow(decay$summary) > 0) {
      out_pdf_decay <- file.path(out_dir, sprintf("LD_decay_chr%s.pdf", chr))
      gg <- ggplot(decay$summary, aes(x = Distance_mid_Mb, y = R2)) +
        geom_point() + geom_line() +
        labs(title = paste0("LD decay — chr ", chr,
                            " (bin ≈ ", format(round(decay$bin_width_bp), big.mark=","), " bp)"),
             x = "Physical distance (Mb)", y = expression(Mean~R^2)) +
        theme_minimal(base_size = 13)
      write.csv(transform(decay$summary, Bin_width_bp = decay$bin_width_bp), file.path(out_dir,paste0("LD_decay_bins_chr",chr,".csv")), row.names=FALSE)
      .efp_ggsave(out_pdf_decay, gg, width = 7.5, height = 5.2, units = "in")
      message("Saved: ", out_pdf_decay)
    } else {
      message("Chr ", chr, ": insufficient pairs for a decay plot after binning.")
    }
  }
}


# =============================================================================
# 16 helpers | Corrected symmetric LD calculations and plotting
# Packages: snpStats, ggplot2, gridExtra, scales. Input schemas: INPUT_FILES.md.
# These functions preserve the supplied latest LD block and its graphics fix.
# =============================================================================
wheat_chromosomes <- as.vector(t(outer(1:7, c("A", "B", "D"), paste0)))
paper_markers <- c("JB1B_662960818", "AX-643813054", "AX-94909517",
                   "AX-95679056", "AX-643821508", "AX-643812441", "PER1")

# Small internal helpers (called automatically by the functions below).
.ld_theme <- function() theme_classic(base_size = 11) +
  theme(axis.ticks = element_line(colour = "black"),
        axis.ticks.length = grid::unit(2, "mm"), legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9))
.save_plot <- function(p, name, out_dir, width = 7, height = 5) {
  # Build lazy .combine(...) arguments BEFORE ggsave opens its output device.
  force(p)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    dest <- file.path(out_dir, paste0(name, ".", ext))
    temp <- tempfile(tmpdir = out_dir, fileext = paste0(".", ext))
    tryCatch({
      ggplot2::ggsave(temp, p, device = ext, width = width, height = height,
             units = "in", dpi = 300, bg = "white", limitsize = FALSE)
      if (!file.exists(temp) || file.info(temp)$size < 100) stop("Empty plot: ", dest)
      if (!file.copy(temp, dest, overwrite = TRUE) || file.info(dest)$size != file.info(temp)$size)
        stop("Could not save: ", dest)
    }, finally = unlink(temp))
  }
}
.combine <- function(plots, ncol = 2, legend = TRUE) {
  # Build on a temporary null device, then restore the device that was active.
  previous_device <- grDevices::dev.cur()
  grDevices::pdf(NULL); device <- grDevices::dev.cur()
  on.exit({
    if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    if (previous_device > 1 && previous_device %in% grDevices::dev.list())
      grDevices::dev.set(previous_device)
  }, add = TRUE)
  key <- NULL
  if (legend) {
    g <- ggplotGrob(plots[[1]] + theme(legend.position = "bottom"))
    for (i in grep("^guide-box", g$layout$name))
      if (!inherits(g$grobs[[i]], "zeroGrob")) { key <- g$grobs[[i]]; break }
  }
  for (i in seq_along(plots)) plots[[i]] <- plots[[i]] +
    labs(tag = letters[i]) + theme(legend.position = "none")
  gridExtra::arrangeGrob(grobs = c(plots, if (!is.null(key)) list(key)), ncol = ncol)
}

# 1. Read, check and align genotypes and map. Return an object used below.
read_ld_data <- function(genotype_file, map_file) {
  read_text <- function(path) {
    if (!file.exists(path)) stop("Input file not found: ", path)
    first <- readLines(path, n = 1, warn = FALSE)
    if (!length(first)) stop("Empty input: ", path)
    read.table(path, header = TRUE, check.names = FALSE, row.names = NULL,
               sep = if (grepl("\t", first, fixed = TRUE)) "\t" else "",
               colClasses = "character", comment.char = "", na.strings = c("NA", "NaN", ".", ""))
  }
  g <- read_text(genotype_file); map <- read_text(map_file)
  if (nrow(g) < 3 || ncol(g) < 3) stop("Need at least 3 samples and 2 SNPs.")
  if (anyNA(g[[1]]) || any(!nzchar(g[[1]])) || anyDuplicated(g[[1]]) || anyDuplicated(names(g)))
    stop("Sample and SNP IDs must be nonempty and unique.")
  raw <- as.matrix(g[-1])
  x <- matrix(suppressWarnings(as.numeric(raw)), nrow = nrow(g),
              dimnames = list(g[[1]], names(g)[-1]))
  if (any(!is.na(raw) & !is.finite(x)) || any(!is.na(x) & !(x %in% 0:2)))
    stop("Genotypes must be 0, 1, 2 or NA; check fractional/missing-value codes.")
  if (!all(c("SNP", "Chromosome", "Position") %in% names(map)) || anyDuplicated(names(map)))
    stop("Map requires unique columns SNP, Chromosome, Position.")
  if (anyNA(map$SNP) || any(!nzchar(map$SNP)) || anyDuplicated(map$SNP)) stop("Invalid/duplicate map SNP IDs.")
  if (length(setdiff(colnames(x), map$SNP))) stop("Some genotype SNPs have no map entry.")
  map <- map[match(colnames(x), map$SNP), c("SNP", "Chromosome", "Position")]
  map$Position <- suppressWarnings(as.numeric(map$Position))
  if (any(!is.finite(map$Position) | map$Position < 0 | map$Position != floor(map$Position)))
    stop("Positions must be nonnegative integer base-pair coordinates.")
  chr <- sub("^CHR", "", toupper(trimws(map$Chromosome)))
  numeric_chr <- !is.na(chr) & grepl("^[0-9]+$", chr)
  chr[numeric_chr] <- wheat_chromosomes[match(chr[numeric_chr], as.character(1:21))]
  if (anyNA(chr) || any(!chr %in% wheat_chromosomes)) stop("Check chromosome codes: expected 1A..7D.")
  map$Chromosome <- chr
  ord <- order(match(chr, wheat_chromosomes), map$Position, map$SNP)
  map <- map[ord, ]; x <- x[, map$SNP, drop = FALSE]; rownames(map) <- NULL
  n <- colSums(!is.na(x))
  usable <- n >= 3 & apply(x, 2, function(v) length(unique(v[!is.na(v)])) >= 2)
  qc <- data.frame(SNP = colnames(x), N_called = n, Retained = usable)
  if (any(!usable)) warning(sum(!usable), " SNPs removed: <3 calls or no genotype variation; see data$qc.")
  list(genotypes = x[, usable, drop = FALSE], map = map[usable, ], qc = qc)
}

# 2. Calculate full symmetric LD ONCE for a chromosome; reuse it for all plots.
# snpStats haplotype r2 is used, with samples in rows and SNPs in columns.
calculate_ld <- function(data, chromosome, max_snps = 1500) {
  map <- data$map[data$map$Chromosome == chromosome, ]
  if (!nrow(map)) stop("No usable SNPs on chromosome ", chromosome)
  if (nrow(map) > max_snps) stop("Increase max_snps explicitly if memory permits full chromosome LD.")
  x <- data$genotypes[, map$SNP, drop = FALSE]; m <- ncol(x)
  observed <- 1 * !is.na(x); n <- crossprod(observed)
  r2 <- if (m == 1) matrix(1, 1, 1) else
    as.matrix(snpStats::ld(methods::as(x, "SnpMatrix"), depth = m-1,
                           stats = "R.squared", symmetric = TRUE))
  xx <- x; xx[is.na(xx)] <- 0
  sums <- crossprod(xx, observed)
  ss <- crossprod(xx*xx, observed) - sums^2/pmax(n, 1)
  r2[n < 3 | ss <= 1e-10 | t(ss) <= 1e-10 | !is.finite(r2)] <- NA
  if (any(r2 < -1e-7 | r2 > 1+1e-7, na.rm = TRUE)) stop("LD outside [0,1].")
  r2[] <- pmax(0, pmin(1, r2)); diag(r2) <- 1  # Self-LD is for display only.
  dimnames(r2) <- list(map$SNP, map$SNP)
  ij <- which(upper.tri(r2), arr.ind = TRUE)   # Count each unordered pair once.
  pairs <- data.frame(SNP1 = map$SNP[ij[,1]], SNP2 = map$SNP[ij[,2]],
                      Distance_bp = abs(map$Position[ij[,2]]-map$Position[ij[,1]]),
                      R2 = r2[ij], N_complete = n[ij])
  list(map = map, r2 = r2, n_complete = n, pairs = pairs)
}

# 3. Lead-marker LD plot and interval table: Fig. 6 (threshold=0.3).
# Also run threshold=0.5 to compare intervals. No automatic threshold relaxation.
plot_marker_ld <- function(ld, marker, threshold = 0.3, out_dir = "LD_outputs") {
  if (length(marker) != 1 || !marker %in% ld$map$SNP) stop("Marker absent from this chromosome: ", marker)
  if (length(threshold) != 1 || !is.finite(threshold) || threshold <= 0 || threshold > 1) stop("Invalid threshold.")
  j <- match(marker, ld$map$SNP)
  d <- data.frame(ld$map, R2 = ld$r2[j, ], N_complete = ld$n_complete[j, ])
  linked <- d$SNP != marker & is.finite(d$R2) & d$R2 >= threshold
  bounds <- range(c(d$Position[j], d$Position[linked]))
  interval <- data.frame(Marker = marker, Chromosome = d$Chromosome[j], Threshold = threshold,
                         Start_bp = bounds[1], End_bp = bounds[2], Width_Mb = diff(bounds)/1e6,
                         Linked_neighbours = sum(linked), Status = if (any(linked)) "LD_envelope" else "collapsed_no_linked_neighbour")
  levels <- c("Lead SNP", paste0("r2 >= ", threshold), paste0("r2 < ", threshold))
  d$Type <- factor(ifelse(d$SNP == marker, levels[1], ifelse(linked, levels[2], levels[3])), levels = levels)
  d <- d[order(d$SNP == marker), ]
  label <- paste0(format(bounds[1], big.mark = ",", scientific = FALSE), " - ",
                  format(bounds[2], big.mark = ",", scientific = FALSE), " bp; ", round(diff(bounds)/1e6, 3),
                  " Mb; ", sum(linked), " linked neighbour(s)")
  p <- ggplot(d, aes(Position/1e6, R2, colour = Type, shape = Type)) +
    geom_vline(xintercept = bounds/1e6, colour = "lightblue3", linetype = "dotted") +
    geom_vline(xintercept = ld$map$Position[j]/1e6, colour = "blue", linetype = "dashed") +
    geom_hline(yintercept = threshold, colour = "#2C9F45", linetype = "dashed") +
    geom_point(size = 1.8, na.rm = TRUE) +
    scale_colour_manual(values = setNames(c("blue", "red", "grey55"), levels), drop = FALSE) +
    scale_shape_manual(values = setNames(c(17, 19, 19), levels), drop = FALSE) +
    scale_x_continuous(breaks = scales::breaks_pretty(5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.03)) +
    labs(title = paste(marker, "|", ld$map$Chromosome[j]),
         subtitle = paste(strwrap(label, 72), collapse = "\n"), x = "Position (Mb)",
         y = expression(LD~(r^2)), colour = NULL, shape = NULL) + .ld_theme()
  name <- paste0("LD_", match(marker, ld$map$SNP), "_", gsub("[^A-Za-z0-9_.-]", "_", marker), "_r2_", threshold)
  .save_plot(p, name, out_dir)
  write.csv(d, file.path(out_dir, paste0(name, "_all_SNPs.csv")), row.names = FALSE)
  write.csv(d[d$SNP != marker & is.finite(d$R2) & d$R2 >= threshold, ],
            file.path(out_dir, paste0(name, "_linked_SNPs.csv")), row.names = FALSE)
  write.csv(interval, file.path(out_dir, paste0(name, "_interval.csv")), row.names = FALSE)
  invisible(list(plot = p, interval = interval))
}

# 4. Chromosome heatmap: Fig. S6. For a local heatmap supply around="PER1".
# SNPs are equally spaced in physical-position order; axis ticks show Mb.
# Grey-blue cells are undefined LD, not zero. Arrows identify supplied markers.
plot_ld_heatmap <- function(ld, markers = character(), around = NULL, window = 10,
                            out_dir = "LD_outputs") {
  if (length(window) != 1 || !is.finite(window) || window < 0 || window != floor(window)) stop("window must be a nonnegative integer.")
  map <- ld$map; r2 <- ld$r2; name <- paste0("LD_heatmap_chr", map$Chromosome[1])
  if (!is.null(around)) {
    j <- match(around, map$SNP)
    if (is.na(j)) stop("Local heatmap marker not found: ", around)
    ii <- seq.int(max(1, j-window), min(nrow(map), j+window))
    map <- map[ii, ]; r2 <- r2[ii, ii, drop = FALSE]; markers <- around
    name <- paste0(name, "_local_", j, "_", gsub("[^A-Za-z0-9_.-]", "_", around))
  }
  m <- nrow(map); ij <- which(upper.tri(r2, diag = TRUE), arr.ind = TRUE)
  d <- data.frame(X = ij[,2], Y = ij[,1], R2 = r2[ij])
  ticks <- unique(round(seq(1, m, length.out = min(5, m))))
  p <- ggplot(d, aes(X, Y, fill = R2)) + geom_tile(width = 1, height = 1) +
    scale_fill_gradientn(colours = c("white", "#FFF0B1", "#F9A13C", "#E33417", "#820000"),
                         limits = c(0, 1), breaks = seq(0, 1, 0.25), na.value = "#B6C1CB") +
    scale_x_continuous(breaks = ticks, labels = round(map$Position[ticks]/1e6, 1)) +
    scale_y_reverse(breaks = ticks, labels = round(map$Position[ticks]/1e6, 1)) + coord_fixed() +
    labs(title = paste("Chromosome", map$Chromosome[1]), x = "Position (Mb; SNP-spaced)",
         y = "Position (Mb; SNP-spaced)", fill = expression(r^2)) + .ld_theme()
  ids <- which(map$SNP %in% markers)
  if (length(ids) && m > 1) {
    a <- data.frame(Label = map$SNP[ids], End = ids, X = 1 + 0.02*(m-1),
                    Y = seq(0.65, 0.9, length.out = length(ids))*(m-1)+1)
    p <- p + geom_segment(data = a, aes(x = X, y = Y-0.025*m, xend = End, yend = End),
                          inherit.aes = FALSE, colour = "red", arrow = grid::arrow(length = grid::unit(1.5, "mm"))) +
      geom_text(data = a, aes(X, Y, label = Label), inherit.aes = FALSE, hjust = 0, vjust = 0, size = 2.5)
  }
  .save_plot(p, name, out_dir, 7, 6)
  invisible(p)
}

# 5. Chromosome-wide LD decay: Fig. S7; excludes self-pairs and undefined LD.
# Adapt bin width until median occupied-bin count reaches min_pairs or the
# width cap is reached. Only bins actually containing >=min_pairs are plotted.
plot_ld_decay <- function(ld, min_pairs = 30, bin_bp = 50000, out_dir = "LD_outputs") {
  if (min_pairs < 1 || bin_bp <= 0) stop("min_pairs and bin_bp must be positive.")
  pairs <- ld$pairs[is.finite(ld$pairs$R2), ]
  b <- data.frame(Bin = integer(), Mean_R2 = numeric(), Pair_count = integer())
  if (nrow(pairs)) {
    cap <- max(bin_bp, 100000, 0.02*max(pairs$Distance_bp))
    repeat {
      bins <- floor(pairs$Distance_bp/bin_bp)
      means <- aggregate(pairs$R2, list(Bin = bins), mean)
      counts <- aggregate(pairs$R2, list(Bin = bins), length)
      b <- data.frame(Bin = means$Bin, Mean_R2 = means$x, Pair_count = counts$x[match(means$Bin, counts$Bin)])
      if (median(b$Pair_count) >= min_pairs || bin_bp >= cap) break
      bin_bp <- min(cap, ceiling(bin_bp*1.5))
    }
  }
  b <- b[order(b$Bin), ]; b$Midpoint_Mb <- (b$Bin+0.5)*bin_bp/1e6
  b$Included <- b$Pair_count >= min_pairs
  d <- b[b$Included, ]; chr <- ld$map$Chromosome[1]
  p <- ggplot(d, aes(Midpoint_Mb, Mean_R2)) + geom_point(size = 1) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    labs(title = paste("Chromosome", chr), subtitle = paste0("Bin width ", signif(bin_bp/1e6, 3),
                                                             " Mb; >=", min_pairs, " pairs/bin"), x = "Physical distance (Mb)", y = expression(Mean~r^2)) + .ld_theme()
  if (nrow(d) > 1) {
    d$Group <- cumsum(c(TRUE, diff(d$Bin) != 1))
    lines <- d[ave(d$Group, d$Group, FUN = length) >= 2, ]
    if (nrow(lines)) p <- p + geom_line(data = lines, aes(group = Group), linewidth = 0.35)
  }
  if (!nrow(d)) p <- p + annotate("text", x = 0, y = 0.5, label = "No bins meet the minimum pair count", size = 3)
  count_plot <- ggplot(b, aes(Midpoint_Mb, Pair_count)) +
    geom_col(width = 0.9*bin_bp/1e6, fill = "skyblue3") + geom_hline(yintercept = min_pairs, linetype = "dashed") +
    labs(title = paste("Pair counts |", chr), x = "Physical distance (Mb)", y = "SNP pairs") + .ld_theme()
  .save_plot(p, paste0("LD_decay_chr", chr), out_dir)
  .save_plot(count_plot, paste0("LD_pair_counts_chr", chr), out_dir)
  b$Bin_width_bp <- rep(bin_bp, nrow(b))
  write.csv(b, file.path(out_dir, paste0("LD_decay_bins_chr", chr, ".csv")), row.names = FALSE)
  invisible(list(plot = p, bins = b))
}

# 6. Optional NEW S5-style diagnostic: Pearson r between consecutive markers.
# This is not the historical GAPIT figure or the haplotype r2 used above.
plot_genome_ld <- function(data, out_dir = "LD_outputs") {
  x <- data$genotypes; map <- data$map; rows <- list()
  for (chr in unique(map$Chromosome)) {
    ids <- which(map$Chromosome == chr)
    if (length(ids) < 2) next
    for (k in seq_len(length(ids)-1)) {
      i <- ids[k]; j <- ids[k+1]; keep <- !is.na(x[,i]) & !is.na(x[,j]); n <- sum(keep)
      r <- if (n >= 3 && sd(x[keep,i]) > 0 && sd(x[keep,j]) > 0) cor(x[keep,i], x[keep,j]) else NA_real_
      rows[[length(rows)+1]] <- data.frame(Chromosome = chr, SNP1 = map$SNP[i], SNP2 = map$SNP[j],
                                           Distance_Mb = (map$Position[j]-map$Position[i])/1e6, R = r, R2 = r^2, N_complete = n)
    }
  }
  if (!length(rows)) stop("No adjacent SNP pairs available.")
  d <- do.call(rbind, rows); d$Index <- seq_len(nrow(d))
  chrs <- unique(d$Chromosome); centres <- sapply(chrs, function(z) mean(range(d$Index[d$Chromosome == z])))
  d$Colour <- factor(match(d$Chromosome, chrs) %% 2)
  genome_axis <- scale_x_continuous(breaks = centres, labels = chrs)
  genome_theme <- theme(axis.text.x = element_text(angle = 90, size = 7), legend.position = "none")
  colours <- scale_colour_manual(values = c("grey60", "#C69518"))
  a <- ggplot(d, aes(Index, R, colour = Colour)) + geom_point(size = 0.6, na.rm = TRUE) +
    genome_axis + colours + labs(x = "Chromosome", y = "Pearson r") + .ld_theme() + genome_theme
  b <- ggplot(d, aes(R)) + geom_histogram(binwidth = 0.1, fill = "grey60", na.rm = TRUE) +
    labs(x = "Pearson r", y = "Pair count") + .ld_theme()
  c <- ggplot(d, aes(Distance_Mb, R)) + geom_point(size = 0.6, na.rm = TRUE) + geom_hline(yintercept = 0, colour = "red") +
    labs(x = "Adjacent-marker distance (Mb)", y = "Pearson r") + .ld_theme()
  e <- ggplot(d, aes(Index, Distance_Mb, colour = Colour)) + geom_point(size = 0.6) +
    genome_axis + colours + labs(x = "Chromosome", y = "Adjacent-marker distance (Mb)") + .ld_theme() + genome_theme
  f <- ggplot(d, aes(Distance_Mb)) + geom_histogram(bins = 40, fill = "grey60") +
    labs(x = "Adjacent-marker distance (Mb)", y = "Pair count") + .ld_theme()
  g <- ggplot(d, aes(Distance_Mb, R2)) + geom_point(size = 0.6, na.rm = TRUE) + coord_cartesian(ylim = c(0,1)) +
    labs(x = "Adjacent-marker distance (Mb)", y = expression(Pearson~r^2)) + .ld_theme()
  ok <- d[is.finite(d$R2), ]
  if (nrow(ok) >= 20 && length(unique(ok$Distance_Mb)) >= 8)
    g <- g + geom_smooth(data = ok, method = "loess", formula = y~x, span = 0.75, se = FALSE, colour = "red")
  p <- .combine(list(a,b,c,e,f,g), 3, legend = FALSE)
  .save_plot(p, "FigS5_adjacent_diagnostics_REANALYSIS", out_dir, 15, 9)
  write.csv(d, file.path(out_dir, "S5_adjacent_Pearson_pairs.csv"), row.names = FALSE)
  invisible(p)
}

# 7. Run everything, analysing each chromosome once and every requested marker.
# All nine markers get individual plots; the seven paper_markers form Fig. 6.
# Results go directly to out_dir (no timestamp folders). A rerun replaces them.
# resume=TRUE redraws from LD_chr*.rds saved by a failed run; use the SAME
# genotype inputs and output folder. It does not recalculate LD.
# New caches also retain the input calls so later resumes can verify them.
run_all_ld <- function(data, markers, out_dir = "LD_outputs", thresholds = c(0.3, 0.5),
                       figure_markers = paper_markers, include_s5 = TRUE, resume = FALSE) {
  markers <- unique(trimws(as.character(markers)))
  if (!length(markers) || anyNA(markers) || any(!nzchar(markers))) stop("Invalid marker list.")
  absent <- setdiff(markers, data$map$SNP)
  if (length(absent)) stop("Requested markers missing/unusable: ", paste(absent, collapse = ", "))
  if (length(setdiff(figure_markers, markers))) stop("All figure_markers must be included in markers.")
  if (!length(figure_markers) || anyDuplicated(figure_markers)) stop("Invalid figure_markers.")
  if (!length(thresholds) || any(!is.finite(thresholds) | thresholds <= 0 | thresholds > 1)) stop("Invalid thresholds.")
  thresholds <- unique(thresholds)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  chr <- data$map$Chromosome[match(markers, data$map$SNP)]
  chromosomes <- wheat_chromosomes[wheat_chromosomes %in% chr]
  if (resume) {
    needed <- file.path(out_dir, paste0("LD_chr", chromosomes, ".rds"))
    if (any(!file.exists(needed))) stop("Saved LD files missing. Run with resume=FALSE.")
    message("Resuming from saved LD. Inputs must be unchanged from the failed run.")
  }
  heatmaps <- decayplots <- intervals <- list()
  plots <- setNames(vector("list", length(thresholds)), as.character(thresholds))
  write.csv(data$qc, file.path(out_dir, "Genotype_QC.csv"), row.names = FALSE)
  write.csv(data$map, file.path(out_dir, "Sorted_map.csv"), row.names = FALSE)
  for (chromosome in chromosomes) {
    message("Chromosome ", chromosome)
    cache_file <- file.path(out_dir, paste0("LD_chr", chromosome, ".rds"))
    if (resume) {
      ld <- readRDS(cache_file)
      expected_map <- data$map[data$map$Chromosome == chromosome, ]
      if (!all(c("map", "r2", "n_complete", "pairs") %in% names(ld)) ||
          !isTRUE(all.equal(ld$map, expected_map, check.attributes = FALSE)))
        stop("Saved LD/map mismatch on ", chromosome, ". Run with resume=FALSE.")
      current_calls <- data$genotypes[, expected_map$SNP, drop = FALSE]
      if (!identical(dim(ld$r2), c(nrow(expected_map), nrow(expected_map))) ||
          !identical(rownames(ld$r2), expected_map$SNP) ||
          !identical(colnames(ld$r2), expected_map$SNP) ||
          !isTRUE(all.equal(ld$n_complete, crossprod(1 * !is.na(current_calls)))))
        stop("Saved LD/sample-count mismatch on ", chromosome, ". Run with resume=FALSE.")
      if (!is.null(ld$input_genotypes) && !identical(ld$input_genotypes, current_calls))
        stop("Saved genotype calls differ on ", chromosome, ". Run with resume=FALSE.")
      if (is.null(ld$input_genotypes))
        message("Using earlier saved LD for ", chromosome,
                "; exact genotype-call verification is unavailable for this cache.")
    } else {
      ld <- calculate_ld(data, chromosome)
      ld$input_genotypes <- data$genotypes[, ld$map$SNP, drop = FALSE]
      saveRDS(ld, cache_file)
    }
    con <- gzfile(file.path(out_dir, paste0("LD_pairs_chr", chromosome, ".csv.gz")), "wt")
    tryCatch(write.csv(ld$pairs, con, row.names = FALSE), finally = close(con))
    heatmaps[[chromosome]] <- plot_ld_heatmap(ld, figure_markers, out_dir = out_dir)
    decayplots[[chromosome]] <- plot_ld_decay(ld, out_dir = out_dir)$plot
    for (marker in markers[chr == chromosome]) {
      plot_ld_heatmap(ld, around = marker, out_dir = out_dir)
      for (threshold in thresholds) {
        r <- plot_marker_ld(ld, marker, threshold, out_dir)
        plots[[as.character(threshold)]][[marker]] <- r$plot
        intervals[[length(intervals)+1]] <- r$interval
      }
    }
  }
  results <- do.call(rbind, intervals)
  write.csv(results, file.path(out_dir, "LD_intervals_all_markers.csv"), row.names = FALSE)
  for (threshold in thresholds)
    .save_plot(.combine(plots[[as.character(threshold)]][figure_markers], 2),
               paste0("Fig6_LD_r2_", threshold), out_dir, 12, 15)
  # Supplementary figure panels use only chromosomes represented by figure_markers.
  figure_chr <- wheat_chromosomes[wheat_chromosomes %in% data$map$Chromosome[match(figure_markers, data$map$SNP)]]
  .save_plot(.combine(heatmaps[figure_chr], 3), "FigS6_LD_heatmaps", out_dir, 15, 10)
  .save_plot(.combine(decayplots[figure_chr], 3, legend = FALSE), "FigS7_LD_decay", out_dir, 15, 9)
  if (include_s5) plot_genome_ld(data, out_dir)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  message("Finished. Outputs: ", normalizePath(out_dir, winslash = "/"))
  invisible(results)
}



# =============================================================================
# 16 | Lead-marker intervals and corrected LD figures
# Packages: snpStats, ggplot2, gridExtra, scales.
# Inputs: numeric genotype/map from 07a; list_of_sigMarkers.csv from 11 (or supplied).
# Outputs: outputs/GWAS/LD_related_plots/; Fig6_LD_r2_0.3/0.5.{pdf,png}, interval/pair/bin CSVs, RDS caches and supplementary alternatives.
# Paper: Fig. 6; S6/S7 alternatives. Optional adjacent-Pearson diagnostic is NOT historical S5.
# Analysis section: 16 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_16 <- function(project_root = ".", resume = FALSE, include_s5_diagnostic = FALSE) {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("snpStats", "ggplot2", "gridExtra", "scales"))
  data <- read_ld_data("outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.Numerical.txt",
                       "outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.map.txt")
  marker_table <- .efp_csv("outputs/GWAS/list_of_sigMarkers.csv", "sigMarker")
  # Both thresholds are exported; neither is relaxed automatically.
  run_all_ld(data, marker_table$sigMarker, out_dir = "outputs/GWAS/LD_related_plots",
             thresholds = c(0.3,0.5), include_s5 = include_s5_diagnostic, resume = resume)
}


# =============================================================================
# 17 | Candidate-gene expression time courses
# Packages: dplyr, tidyr, readr, tibble, ggplot2.
# Inputs: outputs/GWAS/Candidate_genes/Expression_time_course_data.csv; sample, TimeCourse1, TimeCourse2, six gene TPM columns.
# Outputs: outputs/GWAS/Candidate_genes/Expression_time_course_plots11/strongest_candidate_gene_from_each_QTL_expression_corrected.{pdf,png} and _summary.csv.
# Paper: Fig. 7; raw TPM, seven stages; pools TimeCourse2 as in source. SE is undefined for n<2.
# Analysis section: 17 (see RELEASE_NOTES.md for method notes).
# =============================================================================
efp_17 <- function(project_root = ".") {
  project_root <- .efp_root(project_root)
  previous_dir <- getwd(); on.exit(setwd(previous_dir), add = TRUE)
  setwd(project_root)
  .efp_packages(c("dplyr", "tidyr", "readr", "tibble", "ggplot2"))
  # Input and output
  # ============================================================

  in_file <- paste0(
    "outputs/GWAS/Candidate_genes/",
    "Expression_time_course_data.csv"
  )

  out_dir <- paste0(
    "outputs/GWAS/Candidate_genes/",
    "Expression_time_course_plots11"
  )

  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  expr_raw <- .efp_csv(in_file, c("sample","TimeCourse1","TimeCourse2"))
  .efp_unique_ids(expr_raw$sample, "Expression samples")

  # ============================================================
  # Developmental-stage order
  # ============================================================

  time_levels <- c(
    "04 dpa",
    "6-10 dpa",
    "11-15 dpa",
    "16-20 dpa",
    "21-25 dpa",
    "26-30 dpa",
    "35 dpa"
  )

  time_midpoint <- tibble(
    TimeCourse1 = time_levels,
    TimePoint = c(4, 8, 13, 18, 23, 28, 35)
  )

  # ============================================================
  # Strongest candidate gene from each QTL
  #
  # Modify this table if you wish to change any selected gene
  # ============================================================

  candidate_genes <- tribble(
    ~Gene,                 ~QTL,                    ~Function,
    "TraesCS1B02G424500",  "1B: JB1B_662960818",   "XTH/GH16",
    "TraesCS3A02G439700",  "3A: AX-643813054",     "RING/U-box protein",
    "TraesCS5B02G356300",  "5B: AX-94909517",      "UGPase",
    "TraesCS5B02G422600",  "5B: AX-95679056",      "PPR-containing protein",
    "TraesCS6B02G042500",  "6B: PER1",             "Cell-wall peroxidase",
    "TraesCS7A02G349000",  "7A: AX-643812441",     "TBL/DUF231 protein"
  )

  gene_order <- candidate_genes$Gene

  # ============================================================
  # Check that selected genes are present
  # ============================================================

  unknown_stages <- setdiff(na.omit(unique(expr_raw$TimeCourse1)), time_levels)
  if (length(unknown_stages)) stop("Unexpected TimeCourse1 labels: ", paste(unknown_stages, collapse=", "))
  if (anyNA(expr_raw$TimeCourse1)) stop("TimeCourse1 is required for every sample.")
  missing_genes <- setdiff(
    gene_order,
    names(expr_raw)
  )

  if (length(missing_genes) > 0) {
    stop(
      paste0(
        "The following genes were not found in the input file:\n",
        paste(missing_genes, collapse = "\n")
      )
    )
  }

  for (gene in intersect(gene_order, names(expr_raw))) {
    expr_raw[[gene]] <- .efp_number(expr_raw[[gene]], gene)
    if (any(expr_raw[[gene]] < 0, na.rm=TRUE)) stop("Negative TPM in ", gene)
  }

  candidate_genes_present <- candidate_genes %>%
    filter(Gene %in% names(expr_raw))

  if (nrow(candidate_genes_present) == 0) {
    stop("None of the selected genes were found in the input file.")
  }

  # ============================================================
  # Convert selected expression columns to long format
  #
  # TimeCourse2 is retained as metadata but is NOT used to split
  # the developmental expression trajectory.
  # ============================================================

  expr_selected_long <- expr_raw %>%
    select(
      sample,
      TimeCourse1,
      TimeCourse2,
      all_of(candidate_genes_present$Gene)
    ) %>%
    pivot_longer(
      cols = all_of(candidate_genes_present$Gene),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    mutate(
      Expression = .efp_number(Expression, "TPM"),
      TimeCourse1 = factor(
        TimeCourse1,
        levels = time_levels,
        ordered = TRUE
      )
    ) %>%
    left_join(
      time_midpoint,
      by = "TimeCourse1"
    ) %>%
    left_join(
      candidate_genes_present,
      by = "Gene"
    )

  # ============================================================
  # Summarise expression
  #
  # Important:
  # Group only by gene and developmental stage.
  # Do not group by TimeCourse2.
  # ============================================================

  expr_selected_summary <- expr_selected_long %>%
    group_by(
      Gene,
      QTL,
      Function,
      TimeCourse1,
      TimePoint
    ) %>%
    summarise(
      n = sum(!is.na(Expression)),
    
      mean_expr = ifelse(
        n > 0,
        mean(Expression, na.rm = TRUE),
        NA_real_
      ),
    
      sd_expr = ifelse(
        n > 1,
        sd(Expression, na.rm = TRUE),
        NA_real_
      ),
    
      se_expr = ifelse(
        n > 1,
        sd_expr / sqrt(n),
        NA_real_
      ),
    
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(mean_expr - se_expr, 0),
      ymax = mean_expr + se_expr
    ) %>%
    filter(
      !is.na(mean_expr),
      !is.na(TimePoint)
    )

  # ============================================================
  # Set gene and facet order
  # ============================================================

  facet_labels <- candidate_genes_present %>%
    mutate(
      Facet_label = paste0(
        Gene #,
  #      "\n",
  #      QTL,
  #      " | ",
  #      Function
      )
    )

  expr_selected_summary <- expr_selected_summary %>%
    left_join(
      facet_labels %>%
        select(Gene, Facet_label),
      by = "Gene"
    ) %>%
    mutate(
      Gene = factor(
        Gene,
        levels = candidate_genes_present$Gene
      ),
      Facet_label = factor(
        Facet_label,
        levels = facet_labels$Facet_label
      )
    )

  # ============================================================
  # Define colours
  #
  # One colour per candidate gene
  # ============================================================

  gene_colours <- c(
    "TraesCS1B02G424500" = "#0072B2",
    "TraesCS3A02G439700" = "#D55E00",
    "TraesCS5B02G356300" = "#009E73",
    "TraesCS5B02G422600" = "#CC79A7",
    "TraesCS6B02G042500" = "#E69F00",
    "TraesCS7A02G349000" = "#6A3D9A"
  )

  # Retain colours only for genes present
  gene_colours <- gene_colours[
    names(gene_colours) %in% candidate_genes_present$Gene
  ]

  # ============================================================
  # Faceted plot
  #
  # Each facet contains:
  # - one continuous line across seven stages
  # - one ribbon representing ± SE
  # - one gene-specific colour
  # ============================================================

  p_top_candidates <- ggplot(
    expr_selected_summary,
    aes(
      x = TimePoint,
      y = mean_expr,
      colour = Gene,
      fill = Gene,
      group = Gene
    )
  ) +
    geom_ribbon(
      aes(
        ymin = ymin,
        ymax = ymax
      ),
      alpha = 0.22,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1.1,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      size = 2.5,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ Facet_label,
      scales = "free_y",
      ncol = 2
    ) +
    scale_colour_manual(
      values = gene_colours
    ) +
    scale_fill_manual(
      values = gene_colours
    ) +
    scale_x_continuous(
      breaks = time_midpoint$TimePoint,
      labels = time_midpoint$TimeCourse1,
      limits = range(time_midpoint$TimePoint),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(
      title = "Expression profiles of prioritised candidate genes",
      subtitle = "The highest-priority candidate from each gene-containing QTL",
      x = "Developing grain stage",
      y = "Mean expression (TPM)",
      caption = "Line and points represent mean expression; ribbon represents mean ± SE."
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
      axis.text.y = element_text(
        size = 9
      ),
      axis.title = element_text(
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold",
        size = 9,
        lineheight = 1.05
      ),
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40"
      ),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.position = "none"
    )

  print(p_top_candidates)

  # ============================================================
  # Save plot
  # ============================================================

  .efp_ggsave(
    filename = file.path(
      out_dir,
      "strongest_candidate_gene_from_each_QTL_expression_corrected.pdf"
    ),
    plot = p_top_candidates,
    width = 6,
    height = 6
  )

  .efp_ggsave(
    filename = file.path(
      out_dir,
      "strongest_candidate_gene_from_each_QTL_expression_corrected.png"
    ),
    plot = p_top_candidates,
    width = 6,
    height = 6,
    dpi = 600
  )

  # ============================================================
  # Export plotted summary values
  # ============================================================

  write_csv(
    expr_selected_summary,
    file.path(
      out_dir,
      "strongest_candidate_gene_from_each_QTL_expression_corrected_summary.csv"
    )
  )
  invisible(expr_selected_summary)
}

# =============================================================================
# Run selected sections and record the environment
# =============================================================================
efp_run <- function(project_root = ".", sections = c("01", "02", "03", "04")) {
  root <- .efp_root(project_root)
  available <- c(sprintf("%02d", 1:6), "07a", "07b", sprintf("%02d", 8:17))
  if (identical(sections, "all")) sections <- available
  if (anyDuplicated(sections)) {
    message("Repeated section IDs removed; each requested section will run once.")
    sections <- unique(sections)
  }
  if (!length(sections) || any(!sections %in% available))
    stop("Use section IDs: ", paste(available, collapse = ", "))
  folder <- .efp_dir(file.path(root, "outputs/Reproducibility"))
  records <- list()
  on.exit({
    writeLines(capture.output(sessionInfo()), file.path(folder, "sessionInfo.txt"))
    writeLines(capture.output(dput(efp_settings)), file.path(folder, "analysis_settings.R.txt"))
    p <- installed.packages()
    write.csv(p[, intersect(c("Package", "Version", "Built"), colnames(p)), drop=FALSE],
              file.path(folder, "installed_packages.csv"), row.names=FALSE)
    if (length(records)) write.csv(do.call(rbind,records), file.path(folder,"last_run_sections.csv"), row.names=FALSE)
  }, add=TRUE)
  for (id in sections) {
    message("Starting EFP section ", id)
    start <- Sys.time()
    result <- tryCatch({ get(paste0("efp_",id), mode="function")(root); NULL }, error=identity)
    records[[id]] <- data.frame(Section=id, Started=as.character(start),
      Finished=as.character(Sys.time()), Status=if(inherits(result,"error")) "FAILED" else "OK",
      Message=if(inherits(result,"error")) conditionMessage(result) else "")
    if (inherits(result,"error")) stop("Section ",id,": ",conditionMessage(result), call.=FALSE)
  }
  message("Finished. See ",file.path(root,"outputs"))
  invisible(do.call(rbind,records))
}

# =============================================================================
# EXAMPLES -- run the lines you need in R or RStudio
# =============================================================================
# 1. Put this script in your project folder alongside inputs/ and outputs/.
# 2. Load the functions:
# source("EFP_analysis.R")
# 3. Open the repository folder as the R working directory, then use:
# project_root <- "."
# 4. Generate the phenotype figures and summaries:
# efp_run(project_root, c("01", "02", "03", "04"))
#
# Individual sections (each line generates its documented outputs):
# efp_01(project_root)   # Fig. 1a: ordered genotype bar plots
# efp_02(project_root)   # Fig. 1c: histograms
# efp_03(project_root)   # Fig. 1d: correlation plot and statistics
# efp_04(project_root)   # H2 and site CV summaries
# efp_05(project_root)   # Filter single-base HapMap
# efp_06(project_root)   # Filter double-base HapMap
# efp_07a(project_root)  # Convert HapMap; GAPIT genotype diagnostics
# efp_07b(project_root)  # Fig. 2d-e: 100-replicate simulations (long run)
# efp_08(project_root)   # Fig. 2f-g and supplementary PC comparisons (long run)
# efp_09(project_root)   # Final GWAS: MLMM, FarmCPU, BLINK; PCs=0
# efp_10(project_root)   # Fig. 3: Manhattan and QQ plots from saved GWAS
# efp_11(project_root)   # Significant SNPs and correctly labelled allele matrix
# efp_12(project_root)   # Fig. 4: allele effects and ANOVA/Tukey statistics
# efp_13(project_root)   # Figs. 5a/S4: allele-frequency tables and plots
# efp_14(project_root)   # Fig. 5b: favourable-allele count groups
# efp_15(project_root)   # S6/S7: original supplementary layout/binning
# efp_16(project_root)   # Fig. 6 and alternative LD heatmaps/decay plots
# efp_17(project_root)   # Fig. 7: candidate-gene expression
#
# Resume the corrected LD figure block using UNCHANGED saved inputs:
# efp_16(project_root, resume = TRUE)
# Optional NEW adjacent-marker diagnostic; not an exact reproduction of S5:
# efp_16(project_root, resume = TRUE, include_s5_diagnostic = TRUE)
# Reuse complete PC runs only if inputs and GAPIT source have not changed:
# efp_08(project_root, reuse_existing = TRUE)
# To diagnose a GAPIT error, run only the first model/PC combination:
# efp_08(project_root, models = "MLM", pcs = 0)
# This writes to PC_empirical_comparison/diagnostic_subset/ and does not
# overwrite the full paper comparison. A failed run prints its diagnostic log path.
# Entire workflow (requires ALL inputs, including manually prepared tables):
# efp_run(project_root, "all")
#
# Direct LD helper examples (section 16 loads these packages for you):
# .efp_packages(c("snpStats", "ggplot2", "gridExtra", "scales"))
# ld_data <- read_ld_data(
#   file.path(project_root,"outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.Numerical.txt"),
#   file.path(project_root,"outputs/GWAS/Convert_HapMap_format_to_numerical/GAPIT.Genotype.map.txt"))
# ld6B <- calculate_ld(ld_data, "6B")
# plot_marker_ld(ld6B, "PER1", threshold = 0.3, out_dir = file.path(project_root,"outputs/LD_example"))
# plot_ld_heatmap(ld6B, markers = "PER1", out_dir = file.path(project_root,"outputs/LD_example"))
# plot_ld_heatmap(ld6B, around = "PER1", window = 10, out_dir = file.path(project_root,"outputs/LD_example"))
# plot_ld_decay(ld6B, out_dir = file.path(project_root,"outputs/LD_example"))
# plot_genome_ld(ld_data, out_dir = file.path(project_root,"outputs/LD_example"))
#
message("EFP functions loaded. To GENERATE outputs, run efp_run(project_root, c('01','02','03','04')) or a section example above.")
