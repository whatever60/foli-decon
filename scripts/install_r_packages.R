#!/usr/bin/env Rscript

# Configure CRAN mirror.
options(repos = c(CRAN = "https://cloud.r-project.org"))

required_conda_packages <- c(
  "xCell",
  "EPIC",
  "dtangle",
  "BayesPrism",
  "DeconRNASeq",
  "BisqueRNA",
  "DWLS",
  "MuSiC",
  "immunedeconv"
)

missing_conda_packages <- setdiff(required_conda_packages, rownames(installed.packages()))
if (length(missing_conda_packages) > 0) {
  stop(
    paste0(
      "Missing conda-installed R packages: ",
      paste(missing_conda_packages, collapse = ", "),
      ". Recreate/update the conda env from environment.yml before running this script."
    )
  )
}

if (!("remotes" %in% rownames(installed.packages()))) {
  install.packages("remotes")
}

if (!("BiocManager" %in% rownames(installed.packages()))) {
  install.packages("BiocManager")
}

if (!("AnnotationHub" %in% rownames(installed.packages()))) {
  BiocManager::install("AnnotationHub", ask = FALSE, update = FALSE)
}

if (!("xCell2" %in% rownames(installed.packages()))) {
  remotes::install_github(
    "AlmogAngel/xCell2",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}

tracked_packages <- c(required_conda_packages, "xCell2")
installed <- as.data.frame(installed.packages()[, c("Package", "Version")])
print(installed[installed$Package %in% tracked_packages, ])
