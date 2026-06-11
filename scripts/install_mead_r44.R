#!/usr/bin/env Rscript

# Install MEAD packages that are not available from conda-forge/bioconda.
options(repos = c(CRAN = "https://cloud.r-project.org"))

remotes::install_github(
  "mengyin/vashr",
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE
)
remotes::install_github(
  "DongyueXie/MEAD",
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE
)
