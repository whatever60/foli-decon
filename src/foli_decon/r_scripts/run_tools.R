#!/usr/bin/env Rscript

# Read command arguments for tool execution.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript run_tools.R <tool> <payload_json_path>")
}

tool <- args[[1]]
payload_path <- args[[2]]

# Read payload from JSON file.
payload <- jsonlite::fromJSON(payload_path, simplifyVector = FALSE)

# Check whether a payload value is provided as a non-empty string.
is_present <- function(value) {
  !is.null(value) && nzchar(as.character(value))
}

# Convert optional JSON list or scalar into character vector.
to_character_vector <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  as.character(unlist(value))
}

# Read a gene-by-sample matrix from a TSV file.
read_matrix <- function(path) {
  as.matrix(read.table(
    path,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  ))
}

# Read a one-column vector table from TSV.
read_vector <- function(path) {
  table <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  as.vector(table[[1]])
}

# Write a sample-by-celltype matrix to TSV.
write_output <- function(matrix, path) {
  write.table(
    as.data.frame(matrix),
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

# Run xCell through immunedeconv.
run_xcell <- function(payload) {
  options(matrixStats.useNames.NA = "deprecated")
  data("xCell.data", package = "xCell")
  mixture <- read_matrix(payload$mixture_path)
  arrays <- as.logical(payload$arrays)
  expected_cell_types <- to_character_vector(payload$expected_cell_types)

  result <- immunedeconv::deconvolute_xcell(
    gene_expression_matrix = mixture,
    arrays = arrays,
    expected_cell_types = expected_cell_types
  )
  write_output(t(as.matrix(result)), payload$output_path)
}

# Run xCell2 with a pre-trained xCell2 object.
run_xcell2 <- function(payload) {
  library(xCell2)
  source(payload$xcell2_analysis_path)
  mixture <- read_matrix(payload$mixture_path)
  xcell2_object <- readRDS(payload$xcell2_object_path)

  result <- xCell2Analysis(
    mix = mixture,
    xcell2object = xcell2_object,
    minSharedGenes = as.numeric(payload$min_shared_genes),
    rawScores = as.logical(payload$raw_scores),
    spillover = as.logical(payload$spillover),
    spilloverAlpha = as.numeric(payload$spillover_alpha)
  )
  write_output(t(as.matrix(result)), payload$output_path)
}

# Run EPIC with either default or custom signatures.
run_epic <- function(payload) {
  library(EPIC)
  mixture <- read_matrix(payload$mixture_path)
  tumor <- as.logical(payload$tumor)
  scale_mrna <- as.logical(payload$scale_mrna)

  if (is_present(payload$signature_path)) {
    signature <- read_matrix(payload$signature_path)
    signature_genes <- rownames(signature)
    result <- immunedeconv::deconvolute_epic_custom(
      gene_expression_matrix = mixture,
      signature_matrix = signature,
      signature_genes = signature_genes
    )
  } else {
    result <- immunedeconv::deconvolute_epic(
      gene_expression_matrix = mixture,
      tumor = tumor,
      scale_mrna = scale_mrna
    )
  }

  write_output(t(as.matrix(result)), payload$output_path)
}

# Run dtangle with a custom signature matrix.
run_dtangle <- function(payload) {
  mixture <- t(read_matrix(payload$mixture_path))
  signature <- t(read_matrix(payload$signature_path))
  n_markers <- as.numeric(payload$n_markers)
  marker_method <- as.character(payload$marker_method)

  result <- dtangle::dtangle(
    Y = mixture,
    references = signature,
    n_markers = n_markers,
    data_type = "rna-seq",
    marker_method = marker_method
  )
  write_output(as.matrix(result$estimates), payload$output_path)
}

# Run DeconRNASeq with gene-aligned inputs.
run_deconrnaseq <- function(payload) {
  mixture <- as.data.frame(read_matrix(payload$mixture_path), check.names = FALSE)
  signature <- as.data.frame(read_matrix(payload$signature_path), check.names = FALSE)
  use_scale <- as.logical(payload$use_scale)

  if (as.logical(payload$checksig)) {
    message("checksig is currently ignored in the foli-decon DeconRNASeq wrapper.")
  }

  if (use_scale) {
    signature_matrix <- scale(signature)
  } else {
    signature_matrix <- as.matrix(signature)
  }

  out <- matrix(
    NA_real_,
    nrow = ncol(mixture),
    ncol = ncol(signature),
    dimnames = list(colnames(mixture), colnames(signature))
  )

  equality_matrix <- rep(1, ncol(signature_matrix))
  equality_value <- 1
  inequality_matrix <- diag(nrow = ncol(signature_matrix))
  inequality_value <- rep(0, ncol(signature_matrix))

  for (i in seq_len(ncol(mixture))) {
    bulk_sample <- mixture[, i]
    if (use_scale) {
      bulk_sample <- as.vector(scale(bulk_sample))
    } else {
      bulk_sample <- as.vector(bulk_sample)
    }

    solution <- limSolve::lsei(
      A = signature_matrix,
      B = bulk_sample,
      E = equality_matrix,
      F = equality_value,
      G = inequality_matrix,
      H = inequality_value
    )
    out[i, ] <- solution$X
  }

  write_output(out, payload$output_path)
}

# Run DWLS using a provided signature matrix.
run_dwls <- function(payload) {
  bulk <- read_matrix(payload$mixture_path)
  signature <- read_matrix(payload$signature_path)
  method <- as.character(payload$dwls_submethod)

  genes <- intersect(rownames(signature), rownames(bulk))
  bulk <- bulk[genes, , drop = FALSE]
  signature <- signature[genes, , drop = FALSE]

  solutions <- NULL
  for (i in 1:ncol(bulk)) {
    bulk_i <- bulk[, i]

    if (method == "OLS") {
      solution <- DWLS::solveOLS(signature, bulk_i, FALSE)
    } else if (method == "SVR") {
      solution <- DWLS::solveSVR(signature, bulk_i, FALSE)
    } else if (method == "DampenedWLS") {
      solution <- DWLS::solveDampenedWLS(S = signature, B = bulk_i, verbose = FALSE)
    } else {
      stop("dwls_submethod must be one of OLS, SVR, DampenedWLS")
    }

    solutions <- cbind(solutions, solution)
  }

  colnames(solutions) <- colnames(bulk)
  out <- t(solutions)
  write_output(out, payload$output_path)
}

# Run MuSiC using bulk and single-cell count matrices.
run_music <- function(payload) {
  library(Biobase)
  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  batch_ids <- as.character(read_vector(payload$batch_ids_path))
  select_ct <- to_character_vector(payload$select_ct)

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk)
  scrna_eset <- Biobase::ExpressionSet(assayData = scrna)
  Biobase::pData(scrna_eset)$cellType <- cell_types
  Biobase::pData(scrna_eset)$sampleID <- batch_ids

  result <- MuSiC::music_prop(
    bulk.eset = bulk_eset,
    sc.eset = scrna_eset,
    clusters = "cellType",
    samples = "sampleID",
    select.ct = select_ct,
    verbose = FALSE
  )

  write_output(as.matrix(result$Est.prop.weighted), payload$output_path)
}

# Run MuSiC using a full SingleCellExperiment RDS reference.
run_music_sce <- function(payload) {
  library(SingleCellExperiment)
  bulk <- read_matrix(payload$mixture_path)
  sc_sce <- readRDS(payload$scrna_sce_rds_path)
  select_ct <- to_character_vector(payload$select_ct)

  result <- MuSiC::music_prop(
    bulk.mtx = bulk,
    sc.sce = sc_sce,
    clusters = payload$cell_type_column,
    samples = payload$batch_id_column,
    select.ct = select_ct,
    verbose = FALSE
  )

  write_output(as.matrix(result$Est.prop.weighted), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    dir.create(payload$extra_output_dir, showWarnings = FALSE, recursive = TRUE)

    write_df_tsv <- function(x, path) {
      if (is.matrix(x)) {
        x <- as.data.frame(x, check.names = FALSE)
      }
      write.table(
        x,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
      )
    }

    write_named_numeric_tsv <- function(x, path) {
      nm <- names(x)
      if (is.null(nm)) {
        nm <- as.character(seq_along(x))
      }
      df <- data.frame(name = nm, value = as.numeric(x), stringsAsFactors = FALSE)
      write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    write_df_tsv(result[["Est.prop.weighted"]], file.path(payload$extra_output_dir, "Est.prop.weighted.tsv"))
    write_df_tsv(result[["Est.prop.allgene"]], file.path(payload$extra_output_dir, "Est.prop.allgene.tsv"))
    write_df_tsv(result[["Weight.gene"]], file.path(payload$extra_output_dir, "Weight.gene.tsv"))
    write_named_numeric_tsv(result[["r.squared.full"]], file.path(payload$extra_output_dir, "r.squared.full.tsv"))
    write_df_tsv(result[["Var.prop"]], file.path(payload$extra_output_dir, "Var.prop.tsv"))

    manifest <- data.frame(
      key = c("Est.prop.weighted", "Est.prop.allgene", "Weight.gene", "r.squared.full", "Var.prop"),
      type = c("table", "table", "table", "text", "table"),
      file = c("Est.prop.weighted.tsv", "Est.prop.allgene.tsv", "Weight.gene.tsv", "r.squared.full.tsv", "Var.prop.tsv"),
      stringsAsFactors = FALSE
    )
    write.table(
      manifest,
      file = file.path(payload$extra_output_dir, "manifest.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

# Run Bisque reference-based decomposition.
run_bisque <- function(payload) {
  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  batch_ids <- as.character(read_vector(payload$batch_ids_path))

  scrna_eset <- Biobase::ExpressionSet(assayData = scrna)
  Biobase::pData(scrna_eset)$cellType <- cell_types
  Biobase::pData(scrna_eset)$batchId <- batch_ids

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk)

  result <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset = bulk_eset,
    sc.eset = scrna_eset,
    cell.types = "cellType",
    subject.names = "batchId",
    use.overlap = as.logical(payload$use_overlap),
    verbose = FALSE,
    old.cpm = as.logical(payload$old_cpm)
  )

  write_output(t(as.matrix(result$bulk.props)), payload$output_path)
}

# Run BayesPrism with single-cell reference counts.
run_bayesprism <- function(payload) {
  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))

  cell_states <- NULL
  if (is_present(payload$cell_states_path)) {
    cell_states <- as.character(read_vector(payload$cell_states_path))
  }

  update_gibbs <- as.logical(payload$update_gibbs)
  n_cores <- as.integer(payload$n_cores)
  which_theta <- as.character(payload$which_theta)
  state_or_type <- as.character(payload$state_or_type)
  outlier_cut <- as.numeric(payload$outlier_cut)
  outlier_fraction <- as.numeric(payload$outlier_fraction)
  pseudo_min <- as.numeric(payload$pseudo_min)

  if (is.null(cell_states)) {
    cell_states <- cell_types
  }

  prism <- BayesPrism::new.prism(
    reference = t(scrna),
    mixture = t(bulk),
    input.type = "count.matrix",
    cell.type.labels = cell_types,
    cell.state.labels = cell_states,
    key = payload$tum_key,
    outlier.cut = outlier_cut,
    outlier.fraction = outlier_fraction,
    pseudo.min = pseudo_min
  )

  bayes_result <- BayesPrism::run.prism(
    prism = prism,
    n.cores = n_cores,
    update.gibbs = update_gibbs,
    gibbs.control = list(chain.length = 1000, burn.in = 500, thinning = 2),
    opt.control = list(trace = 0, maxit = 100000)
  )

  if (update_gibbs) {
    theta <- BayesPrism::get.fraction(
      bp = bayes_result,
      which.theta = which_theta,
      state.or.type = state_or_type
    )
  } else {
    theta <- bayes_result@posterior.initial.cellType@theta
  }

  write_output(as.matrix(theta), payload$output_path)
}

# Dispatch tool selection.
if (tool == "xcell") {
  run_xcell(payload)
} else if (tool == "xcell2") {
  run_xcell2(payload)
} else if (tool == "epic") {
  run_epic(payload)
} else if (tool == "dtangle") {
  run_dtangle(payload)
} else if (tool == "deconrnaseq") {
  run_deconrnaseq(payload)
} else if (tool == "dwls") {
  run_dwls(payload)
} else if (tool == "music") {
  run_music(payload)
} else if (tool == "music_sce") {
  run_music_sce(payload)
} else if (tool == "bisque") {
  run_bisque(payload)
} else if (tool == "bayesprism") {
  run_bayesprism(payload)
} else {
  stop(paste0("Unsupported R tool: ", tool))
}
