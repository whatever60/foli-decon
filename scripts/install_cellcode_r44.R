#!/usr/bin/env Rscript

# Install CellCODE after conda has provided its compiled/Bioconductor deps.
options(repos = c(CRAN = "https://cloud.r-project.org"))

remotes::install_github(
  "mchikina/CellCODE",
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE
)
