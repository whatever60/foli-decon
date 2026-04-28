#!/usr/bin/env Rscript

required_packages <- c(
  "Biobase",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "TOAST",
  "Matrix",
  "MCMCpack",
  "nnls",
  "ggplot2",
  "remotes"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install missing conda-managed packages first: ", paste(missing_packages, collapse = ", "))
}

remotes::install_github(
  "xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c",
  upgrade = "never",
  dependencies = FALSE,
  build_vignettes = FALSE
)

library(MuSiC)
if (as.character(packageVersion("MuSiC")) != "1.0.0") {
  stop("Expected MuSiC 1.0.0, found ", as.character(packageVersion("MuSiC")))
}
