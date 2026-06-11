required_packages <- c(
  "Biobase",
  "dplyr",
  "tidyr",
  "Matrix",
  "caret",
  "NMF",
  "pheatmap",
  "abind",
  "pbapply",
  "Rcpp",
  "RcppArmadillo"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install missing conda-managed packages first: ", paste(missing_packages, collapse = ", "))
}

remotes::install_github(
  "humengying0907/InstaPrism",
  dependencies = FALSE,
  upgrade = "never"
)
