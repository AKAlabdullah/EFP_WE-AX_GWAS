# Install dependencies once; sourcing this file only defines the installer.
# Run: source("install_packages.R"); efp_install_packages(project_root = ".")
# GAPIT is loaded from vendor/gapit_functions.txt by the analysis script.
# An existing GAPIT source is preserved and NEVER overwritten.
# GAPIT 4.1 dependencies include its optional interactive and parallel routines.

efp_cran_packages <- c("dplyr", "tidyr", "readr", "tibble", "ggplot2", "tidytext",
  "lme4", "rMVP", "gridExtra", "multcompView", "scales",
  "gplots", "genetics", "ape", "bigmemory", "EMMREML", "scatterplot3d",
  "MASS", "RcppArmadillo", "snowfall", "plotly", "htmltools", "magrittr")
efp_bioc_packages <- c("snpStats", "multtest")

efp_package_status <- function() {
  pkgs <- c(efp_cran_packages, efp_bioc_packages)
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  data.frame(Package=pkgs, Installed=ok,
    Version=vapply(pkgs, function(p) if(requireNamespace(p,quietly=TRUE))
      as.character(utils::packageVersion(p)) else NA_character_, character(1)), row.names=NULL)
}

efp_install_packages <- function(project_root = ".", download_gapit = FALSE) {
  if (!dir.exists(project_root)) stop("Project folder not found: ",project_root)
  project_root <- normalizePath(project_root, winslash="/", mustWork=TRUE)
  repo <- "https://cloud.r-project.org"
  need <- efp_cran_packages[!vapply(efp_cran_packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(need)) install.packages(need, repos=repo)
  if (requireNamespace("ggplot2",quietly=TRUE) && utils::packageVersion("ggplot2") < "3.4.0")
    install.packages("ggplot2",repos=repo)
  need <- efp_bioc_packages[!vapply(efp_bioc_packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(need)) {
    if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos=repo)
    BiocManager::install(need, ask=FALSE, update=FALSE)
  }
  gapit_file <- file.path(project_root,"vendor","gapit_functions.txt")
  dir.create(dirname(gapit_file),recursive=TRUE,showWarnings=FALSE)
  if (!file.exists(gapit_file) && download_gapit) {
    # Explicitly requested current-source download, NOT a historical GAPIT pin.
    # Place your study's archived source here instead for exact reproduction.
    url <- "https://zzlab.net/GAPIT/gapit_functions.txt"
    tmp <- tempfile(tmpdir=dirname(gapit_file), fileext=".txt")
    on.exit(unlink(tmp),add=TRUE)
    download.file(url,destfile=tmp,mode="wb")
    if (!file.exists(tmp) || file.info(tmp)$size == 0) stop("GAPIT download is empty.")
    invisible(parse(file=tmp))
    if (!file.copy(tmp,gapit_file,overwrite=FALSE)) stop("Could not preserve downloaded GAPIT source.")
    write.csv(data.frame(URL=url, Downloaded_UTC=format(Sys.time(),tz="UTC",usetz=TRUE),
      MD5=unname(tools::md5sum(gapit_file)), Provenance="New download; historical study version unverified"),
      file.path(dirname(gapit_file),"GAPIT_download_provenance.csv"),row.names=FALSE)
    message("Saved current GAPIT source. This is NOT verified as the version used in the paper.")
  }
  if (file.exists(gapit_file)) {
    Sys.setenv(EFP_GAPIT_SOURCE=normalizePath(gapit_file,winslash="/",mustWork=TRUE))
    message("GAPIT source preserved: ",gapit_file)
  } else {
    message("R dependencies installed. For GWAS, copy the study's GAPIT source to vendor/gapit_functions.txt.\n",
      "For a deliberate rerun with the current source, use download_gapit=TRUE.")
  }
  status <- efp_package_status(); print(status)
  if (any(!status$Installed)) stop("Installation incomplete: ",paste(status$Package[!status$Installed],collapse=", "))
  if (utils::packageVersion("ggplot2") < "3.4.0") stop("ggplot2 still needs updating to >=3.4.0.")
  invisible(status)
}

# EXAMPLES:
# source("install_packages.R")
# efp_install_packages(project_root = ".")
# efp_package_status()
# Current GAPIT download only if no saved copy exists (changes version provenance):
# efp_install_packages(project_root = ".", download_gapit = TRUE)
