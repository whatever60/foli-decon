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

# Convert optional JSON scalar into a numeric value.
to_optional_numeric <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  as.numeric(value)
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

# Read an ABIS signature in either native source format or wrapper TSV format.
read_abis_signature <- function(path) {
  header_fields <- strsplit(readLines(path, n = 1), "\t", fixed = TRUE)[[1]]
  table <- read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  if (ncol(table) == length(header_fields)) {
    return(as.matrix(table))
  }
  rownames(table) <- as.character(table[[1]])
  table <- table[, -1, drop = FALSE]
  as.matrix(table)
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

# Write cell-type-first data.frame as sample-by-celltype TSV.
write_cell_type_first_output <- function(result, path) {
  if (is.data.frame(result) && identical(colnames(result)[[1]], "cell_type")) {
    cell_types <- as.character(result[[1]])
    values <- as.matrix(result[, -1, drop = FALSE])
    mode(values) <- "numeric"
    out <- t(values)
    colnames(out) <- cell_types
    write_output(out, path)
  } else {
    write_output(t(as.matrix(result)), path)
  }
}

# Write a matrix/data.frame to TSV with row names.
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

# Write a GCT v1.2 file for native tools that require GCT input.
write_gct <- function(matrix, path) {
  table <- data.frame(
    NAME = rownames(matrix),
    Description = rownames(matrix),
    as.data.frame(matrix, check.names = FALSE),
    check.names = FALSE
  )
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("#1.2", con = con)
  writeLines(paste(nrow(matrix), ncol(matrix), sep = "\t"), con = con)
  write.table(
    table,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# Read a GCT v1.2 score table and return sample-by-score output.
read_gct_scores <- function(path) {
  table <- read.table(
    path,
    sep = "\t",
    skip = 2,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  rownames(table) <- as.character(table$NAME)
  table$NAME <- NULL
  table$Description <- NULL
  result <- t(as.matrix(table))
  mode(result) <- "numeric"
  result
}

# Write a named vector/list to a simple two-column TSV.
write_named_values_tsv <- function(x, path) {
  values <- unlist(x)
  nm <- names(values)
  if (is.null(nm)) {
    nm <- as.character(seq_along(values))
  }
  df <- data.frame(name = nm, value = as.character(values), stringsAsFactors = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Persist selected diagnostics beside the main matrix output.
write_result_extras <- function(result, extra_output_dir, keys) {
  dir.create(extra_output_dir, showWarnings = FALSE, recursive = TRUE)

  keys <- intersect(keys, names(result))
  for (key in keys) {
    value <- result[[key]]
    path <- file.path(extra_output_dir, paste0(key, ".tsv"))
    if (is.matrix(value) || is.data.frame(value)) {
      write_df_tsv(value, path)
    } else {
      write_named_values_tsv(value, path)
    }
  }

  manifest <- data.frame(
    key = keys,
    file = paste0(keys, ".tsv"),
    stringsAsFactors = FALSE
  )
  write.table(
    manifest,
    file = file.path(extra_output_dir, "manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Persist optional MuSiC2 diagnostics beside the main matrix output.
write_music2_extras <- function(result, extra_output_dir) {
  write_result_extras(
    result = result,
    extra_output_dir = extra_output_dir,
    keys = c("Est.prop", "convergence", "n.iter", "DE.genes", "id.not.converge")
  )
}

# Run xCell through its native R API.
run_xcell <- function(payload) {
  options(matrixStats.useNames.NA = "deprecated")
  data("xCell.data", package = "xCell")
  mixture <- read_matrix(payload$mixture_path)
  arrays <- as.logical(payload$arrays)
  expected_cell_types <- to_character_vector(payload$expected_cell_types)

  result <- xCell::xCellAnalysis(
    expr = mixture,
    rnaseq = !arrays,
    cell.types.use = expected_cell_types
  )
  write_output(t(as.matrix(result$Scores)), payload$output_path)
}

# Run quanTIseq through its native R API.
run_quantiseq <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  result <- quantiseqr::run_quantiseq(
    expression_data = mixture,
    signature_matrix = "TIL10",
    is_arraydata = as.logical(payload$arrays),
    is_tumordata = as.logical(payload$tumor),
    scale_mRNA = as.logical(payload$scale_mrna)
  )
  result <- as.data.frame(result, check.names = FALSE)
  if ("Sample" %in% colnames(result)) {
    rownames(result) <- as.character(result$Sample)
    result$Sample <- NULL
  }
  write_output(as.matrix(result), payload$output_path)
}

# Run MCP-counter through its native API.
run_mcp_counter <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  mcp_args <- list(
    expression = mixture,
    featuresType = as.character(payload$feature_types)
  )
  if (is_present(payload$probesets_path)) {
    mcp_args$probesets <- read.table(
      payload$probesets_path,
      sep = "\t",
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
  }
  if (is_present(payload$genes_path)) {
    mcp_args$genes <- read.table(
      payload$genes_path,
      sep = "\t",
      stringsAsFactors = FALSE,
      header = TRUE,
      colClasses = "character",
      check.names = FALSE
    )
  }
  result <- do.call(MCPcounter::MCPcounter.estimate, mcp_args)
  write_output(t(as.matrix(result)), payload$output_path)
}

# Run PSEA marker reference-signal scoring.
run_psea <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  markers <- read.table(
    payload$marker_sets_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  sample_subset <- NULL
  if (!is.null(payload$sample_subset)) {
    sample_subset_values <- unlist(payload$sample_subset)
    if (length(sample_subset_values) > 0) {
      sample_subset_values <- as.character(sample_subset_values)
      if (all(grepl("^[0-9]+$", sample_subset_values))) {
        sample_subset <- as.integer(sample_subset_values)
      } else {
        sample_subset <- match(sample_subset_values, colnames(mixture))
        if (any(is.na(sample_subset))) {
          stop("sample_subset contains sample names not found in mixture columns")
        }
      }
    }
  }

  cell_types <- unique(markers$cell_type)
  out <- matrix(
    NA_real_,
    nrow = ncol(mixture),
    ncol = length(cell_types),
    dimnames = list(colnames(mixture), cell_types)
  )

  for (cell_type in cell_types) {
    marker_rows <- markers[markers$cell_type == cell_type, , drop = FALSE]
    grouped_markers <- split(marker_rows$gene, marker_rows$marker_group)
    out[, cell_type] <- PSEA::marker(
      expr = mixture,
      id = grouped_markers,
      sampleSubset = sample_subset,
      targetMean = as.numeric(payload$target_mean)
    )
  }

  write_output(out, payload$output_path)
}

# Run CellCODE surrogate proportion variable scoring.
run_cellcode <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  signature <- read_matrix(payload$signature_path)
  groups <- as.character(read_vector(payload$groups_path))

  max_markers <- NULL
  if (!is.null(payload$max_markers)) {
    max_markers <- as.integer(payload$max_markers)
  }

  data_tag <- CellCODE::tagData(
    data = signature,
    cutoff = as.numeric(payload$tag_cutoff),
    max = max_markers,
    ref = mixture,
    ref.mean = as.logical(payload$ref_mean)
  )
  result <- CellCODE::getAllSPVs(
    data = mixture,
    grp = groups,
    dataTag = data_tag,
    method = as.character(payload$method),
    plot = FALSE,
    mix.par = as.numeric(payload$mix_par)
  )
  rownames(result) <- colnames(mixture)

  write_output(as.matrix(result), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    dir.create(payload$extra_output_dir, showWarnings = FALSE, recursive = TRUE)
    write_df_tsv(data_tag, file.path(payload$extra_output_dir, "data_tag.tsv"))
    write_df_tsv(result, file.path(payload$extra_output_dir, "spv.tsv"))
    manifest <- data.frame(
      key = c("data_tag", "spv"),
      file = c("data_tag.tsv", "spv.tsv"),
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

# Run ABIS using the source logic from the local Shiny app implementation.
run_abis <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  technology <- as.character(payload$technology)
  signature <- read_abis_signature(payload$signature_path)
  genes <- intersect(rownames(mixture), rownames(signature))
  if (length(genes) == 0) {
    stop("ABIS found zero shared genes between mixture and signature")
  }

  if (technology == "rnaseq") {
    mixture_for_fit <- mixture[genes, , drop = FALSE]
  } else if (technology == "microarray") {
    target <- read.table(payload$target_path, header = FALSE)
    mixture_for_fit <- preprocessCore::normalize.quantiles.use.target(
      as.matrix(mixture[genes, , drop = FALSE]),
      target[, 1],
      copy = TRUE,
      subset = NULL
    )
    colnames(mixture_for_fit) <- colnames(mixture)
    rownames(mixture_for_fit) <- genes
  } else {
    stop("ABIS technology must be 'rnaseq' or 'microarray'")
  }

  fit_matrix <- as.matrix(signature[genes, , drop = FALSE])
  coefficients <- sapply(
    seq_len(ncol(mixture_for_fit)),
    function(sample_index) {
      coef(MASS::rlm(
        x = fit_matrix,
        y = mixture_for_fit[, sample_index],
        maxit = as.integer(payload$maxit)
      ))
    }
  )
  if (is.null(dim(coefficients))) {
    coefficients <- matrix(coefficients, ncol = 1)
  }
  rownames(coefficients) <- colnames(signature)
  colnames(coefficients) <- colnames(mixture_for_fit)
  write_output(t(signif(coefficients * 100, 3)), payload$output_path)
}

# Run ESTIMATE through the original R-Forge package API.
run_estimate <- function(payload) {
  library(estimate)
  mixture <- read_matrix(payload$mixture_path)
  input_gct <- tempfile(pattern = "foli_estimate_input_", fileext = ".gct")
  output_gct <- tempfile(pattern = "foli_estimate_output_", fileext = ".gct")
  write_gct(mixture, input_gct)
  estimate::estimateScore(
    input.ds = input_gct,
    output.ds = output_gct,
    platform = as.character(payload$platform)
  )
  write_output(read_gct_scores(output_gct), payload$output_path)
}

# Run Consensus-TME through its native R API.
run_consensus_tme <- function(payload) {
  mixture <- read_matrix(payload$mixture_path)
  cancer_type <- NULL
  if (is_present(payload$cancer_type)) {
    cancer_type <- as.character(payload$cancer_type)
  }
  result <- ConsensusTME::consensusTMEAnalysis(
    bulkExp = mixture,
    cancerType = cancer_type,
    statMethod = as.character(payload$stat_method),
    singScoreDisp = as.logical(payload$sing_score_disp),
    immuneScore = as.logical(payload$immune_score),
    excludeCells = to_character_vector(payload$exclude_cells),
    parallel.sz = as.integer(payload$parallel_size)
  )
  write_output(t(as.matrix(result$Scores)), payload$output_path)
}

# TIMER local native support is deferred.
run_timer <- function(payload) {
  stop("TIMER local native support is deferred; TIMER2/TIMER3 web wrappers are not used.")
}

# Run CDSeq complete deconvolution with optional reference labels.
run_cdseq <- function(payload) {
  bulk <- read_matrix(payload$mixture_path)
  cell_type_number <- as.integer(unlist(payload$cell_type_number))

  beta <- NULL
  if (!is.null(payload$beta)) {
    beta <- as.numeric(unlist(payload$beta))
  }

  gene_subset_size <- NULL
  if (!is.null(payload$gene_subset_size)) {
    gene_subset_size <- as.integer(payload$gene_subset_size)
  }

  cpu_number <- NULL
  if (!is.null(payload$cpu_number)) {
    cpu_number <- as.integer(payload$cpu_number)
  }

  gene_length <- NULL
  if (is_present(payload$gene_lengths_path)) {
    gene_length <- as.numeric(read_vector(payload$gene_lengths_path))
  }

  reference_gep <- NULL
  if (is_present(payload$reference_gep_path)) {
    reference_gep <- read_matrix(payload$reference_gep_path)
  }

  cdseq_args <- list(
    bulk_data = bulk,
    beta = beta,
    alpha = as.numeric(payload$alpha),
    cell_type_number = cell_type_number,
    mcmc_iterations = as.integer(payload$mcmc_iterations),
    dilution_factor = as.numeric(payload$dilution_factor),
    block_number = as.integer(payload$block_number)
  )
  if (!is.null(gene_subset_size)) {
    cdseq_args$gene_subset_size <- gene_subset_size
  }
  if (!is.null(cpu_number)) {
    cdseq_args$cpu_number <- cpu_number
  }
  if (!is.null(gene_length)) {
    cdseq_args$gene_length <- gene_length
  }
  if (!is.null(reference_gep)) {
    cdseq_args$reference_gep <- reference_gep
  }

  result <- do.call(CDSeq::CDSeq, cdseq_args)

  write_output(t(as.matrix(result$estProp)), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    write_result_extras(
      result = result,
      extra_output_dir = payload$extra_output_dir,
      keys = c("estProp", "estGEP", "cell_type_assignment", "cellTypeAssignSplit", "lgpst", "estT")
    )
  }
}

# Run xCell2 with a pre-trained xCell2 object.
run_xcell2 <- function(payload) {
  library(xCell2)
  library(BiocParallel)
  source(payload$xcell2_analysis_path)
  mixture <- read_matrix(payload$mixture_path)
  xcell2_object <- readRDS(payload$xcell2_object_path)
  xcell2_workers <- as.integer(payload$xcell2_workers)
  if (xcell2_workers > 1) {
    bp <- BiocParallel::MulticoreParam(workers = xcell2_workers)
  } else {
    bp <- BiocParallel::SerialParam()
  }

  result <- xCell2Analysis(
    mix = mixture,
    xcell2object = xcell2_object,
    minSharedGenes = as.numeric(payload$min_shared_genes),
    rawScores = as.logical(payload$raw_scores),
    spillover = as.logical(payload$spillover),
    spilloverAlpha = as.numeric(payload$spillover_alpha),
    BPPARAM = bp
  )
  write_cell_type_first_output(result, payload$output_path)
}

# Train a custom xCell2 reference object.
run_xcell2_train <- function(payload) {
  library(xCell2)
  library(BiocParallel)

  ref <- read_matrix(payload$scrna_path)
  labels <- read.table(
    payload$labels_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  xcell2_workers <- as.integer(payload$xcell2_workers)
  if (xcell2_workers > 1) {
    bp <- BiocParallel::MulticoreParam(workers = xcell2_workers)
  } else {
    bp <- BiocParallel::SerialParam()
  }

  lineage_file <- NULL
  if (is_present(payload$lineage_file)) {
    lineage_file <- as.character(payload$lineage_file)
  }

  set.seed(as.integer(payload$seed))
  ref_obj <- xCell2::xCell2Train(
    ref = ref,
    labels = labels,
    refType = as.character(payload$ref_type),
    lineageFile = lineage_file,
    BPPARAM = bp,
    useOntology = as.logical(payload$use_ontology),
    returnSignatures = as.logical(payload$return_signatures),
    returnAnalysis = FALSE,
    useSpillover = as.logical(payload$use_spillover),
    spilloverAlpha = as.numeric(payload$spillover_alpha),
    minPbCells = as.integer(payload$min_pb_cells),
    minPbSamples = as.integer(payload$min_pb_samples),
    minScGenes = as.integer(payload$min_sc_genes)
  )
  saveRDS(ref_obj, file = payload$xcell2_object_path)

  status <- matrix(
    c("ok", payload$xcell2_object_path),
    nrow = 1,
    dimnames = list("xcell2_train", c("status", "xcell2_object_path"))
  )
  write_output(status, payload$output_path)
}

# Run EPIC with either default or custom signatures.
run_epic <- function(payload) {
  library(EPIC)
  mixture <- read_matrix(payload$mixture_path)
  tumor <- as.logical(payload$tumor)
  scale_mrna <- as.logical(payload$scale_mrna)

  if (is_present(payload$signature_path)) {
    signature <- read_matrix(payload$signature_path)
    reference <- list(
      refProfiles = signature,
      sigGenes = rownames(signature)
    )
    result <- EPIC::EPIC(
      bulk = mixture,
      reference = reference,
      scaleExprs = scale_mrna,
      withOtherCells = tumor
    )
  } else {
    result <- EPIC::EPIC(
      bulk = mixture,
      scaleExprs = scale_mrna,
      withOtherCells = tumor
    )
  }

  write_output(as.matrix(result$cellFractions), payload$output_path)
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

    write_df_tsv(result[["Est.prop.weighted"]], file.path(payload$extra_output_dir, "Est.prop.weighted.tsv"))
    write_df_tsv(result[["Est.prop.allgene"]], file.path(payload$extra_output_dir, "Est.prop.allgene.tsv"))
    write_df_tsv(result[["Weight.gene"]], file.path(payload$extra_output_dir, "Weight.gene.tsv"))
    write_named_values_tsv(result[["r.squared.full"]], file.path(payload$extra_output_dir, "r.squared.full.tsv"))
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

# Run MuSiC2 once the R-side reference object is ready.
run_music2_common <- function(payload, bulk_control, bulk_case, sc_sce, clusters, samples) {
  select_ct <- to_character_vector(payload$select_ct)

  result <- MuSiC::music2_prop_t_statistics(
    bulk.control.mtx = bulk_control,
    bulk.case.mtx = bulk_case,
    sc.sce = sc_sce,
    clusters = clusters,
    samples = samples,
    select.ct = select_ct,
    expr_low = as.numeric(payload$expr_low),
    prop_r = as.numeric(payload$prop_r),
    eps_c = as.numeric(payload$eps_c),
    eps_r = as.numeric(payload$eps_r),
    n_resample = as.integer(payload$n_resample),
    sample_prop = as.numeric(payload$sample_prop),
    cutoff_expr = as.numeric(payload$cutoff_expr),
    cutoff_fc = as.numeric(payload$cutoff_fc),
    cutoff_c = as.numeric(payload$cutoff_c),
    cutoff_r = as.numeric(payload$cutoff_r),
    maxiter = as.integer(payload$maxiter),
    ct.cov = as.logical(payload$ct_cov),
    centered = as.logical(payload$centered),
    normalize = as.logical(payload$normalize)
  )

  write_output(as.matrix(result$Est.prop), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    write_music2_extras(result, payload$extra_output_dir)
  }
}

# Run MuSiC2 using bulk and single-cell count matrices.
run_music2 <- function(payload) {
  library(SingleCellExperiment)
  library(SummarizedExperiment)

  bulk_control <- read_matrix(payload$control_mixture_path)
  bulk_case <- read_matrix(payload$case_mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  batch_ids <- as.character(read_vector(payload$batch_ids_path))

  sc_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = scrna))
  SummarizedExperiment::colData(sc_sce)$cellType <- cell_types
  SummarizedExperiment::colData(sc_sce)$sampleID <- batch_ids

  run_music2_common(
    payload = payload,
    bulk_control = bulk_control,
    bulk_case = bulk_case,
    sc_sce = sc_sce,
    clusters = "cellType",
    samples = "sampleID"
  )
}

# Run MuSiC2 using a full SingleCellExperiment RDS reference.
run_music2_sce <- function(payload) {
  library(SingleCellExperiment)

  bulk_control <- read_matrix(payload$control_mixture_path)
  bulk_case <- read_matrix(payload$case_mixture_path)
  sc_sce <- readRDS(payload$scrna_sce_rds_path)

  run_music2_common(
    payload = payload,
    bulk_control = bulk_control,
    bulk_case = bulk_case,
    sc_sce = sc_sce,
    clusters = payload$cell_type_column,
    samples = payload$batch_id_column
  )
}

# Run SCDC using bulk and single-cell count matrices.
run_scdc <- function(payload) {
  library(Biobase)

  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  batch_ids <- as.character(read_vector(payload$batch_ids_path))
  ct_sub <- to_character_vector(payload$ct_sub)
  if (is.null(ct_sub)) {
    ct_sub <- unique(cell_types)
  }

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk)
  scrna_eset <- Biobase::ExpressionSet(assayData = scrna)
  Biobase::pData(scrna_eset)$cellType <- cell_types
  Biobase::pData(scrna_eset)$sampleID <- batch_ids

  result <- SCDC::SCDC_prop(
    bulk.eset = bulk_eset,
    sc.eset = scrna_eset,
    ct.varname = "cellType",
    sample = "sampleID",
    ct.sub = ct_sub,
    iter.max = as.integer(payload$iter_max),
    nu = as.numeric(payload$nu),
    epsilon = as.numeric(payload$epsilon),
    weight.basis = as.logical(payload$weight_basis),
    Transform_bisque = as.logical(payload$transform_bisque)
  )

  write_output(as.matrix(result$prop.est.mvw), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    write_result_extras(
      result = result,
      extra_output_dir = payload$extra_output_dir,
      keys = c("prop.est.mvw", "basis.mvw", "yhat", "yeval", "peval")
    )
  }
}

# Run MEAD using bulk and single-cell count matrices.
run_mead <- function(payload) {
  library(SingleCellExperiment)
  library(SummarizedExperiment)

  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  batch_ids <- as.character(read_vector(payload$batch_ids_path))

  sc_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = scrna))
  SummarizedExperiment::colData(sc_sce)$cell_type <- cell_types
  SummarizedExperiment::colData(sc_sce)$individual <- batch_ids

  groups <- NULL
  if (is_present(payload$groups_path)) {
    groups <- as.character(read_vector(payload$groups_path))
  }

  preprocessing_control <- list(
    marker_gene = to_character_vector(payload$marker_genes),
    gene_thresh = as.numeric(payload$gene_thresh),
    max_count_quantile_celltype = to_optional_numeric(payload$max_count_quantile_celltype),
    max_count_quantile_indi = to_optional_numeric(payload$max_count_quantile_indi),
    filter.gene = as.logical(payload$filter_gene)
  )
  estimation_control <- list(
    hc.type = as.character(payload$hc_type),
    centeringXY = as.logical(payload$centering_xy),
    nfold = as.integer(payload$nfold),
    groups = groups,
    calc_var = as.logical(payload$calc_var),
    use.QP = as.logical(payload$use_qp)
  )

  result <- MEAD::MEAD(
    bulk = bulk,
    ref = sc_sce,
    cell_types = to_character_vector(payload$select_ct),
    preprocessing_control = preprocessing_control,
    estimation_control = estimation_control
  )

  write_output(t(as.matrix(result$p_hat)), payload$output_path)
  if (as.logical(payload$calc_var)) {
    write_output(t(as.matrix(result$p_hat_se)), payload$uncertainty_path)
  }

  if (is_present(payload$extra_output_dir)) {
    dir.create(payload$extra_output_dir, showWarnings = FALSE, recursive = TRUE)
    write_df_tsv(result$p_hat, file.path(payload$extra_output_dir, "p_hat_celltype_by_sample.tsv"))
    if (as.logical(payload$calc_var)) {
      write_df_tsv(result$p_hat_se, file.path(payload$extra_output_dir, "p_hat_se_celltype_by_sample.tsv"))
      ci <- MEAD::get_ci(result)
      write_df_tsv(ci$ci_l, file.path(payload$extra_output_dir, "ci_l.tsv"))
      write_df_tsv(ci$ci_r, file.path(payload$extra_output_dir, "ci_r.tsv"))
    }
    if ("p_group_diff_hat" %in% names(result)) {
      write_df_tsv(result$p_group_diff_hat, file.path(payload$extra_output_dir, "p_group_diff_hat.tsv"))
    }
    manifest <- data.frame(
      key = c("p_hat", "p_hat_se", "ci_l", "ci_r", "p_group_diff_hat"),
      file = c(
        "p_hat_celltype_by_sample.tsv",
        "p_hat_se_celltype_by_sample.tsv",
        "ci_l.tsv",
        "ci_r.tsv",
        "p_group_diff_hat.tsv"
      ),
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

# Run InstaPrism using a reference prepared from single-cell data.
run_instaprism <- function(payload) {
  bulk <- read_matrix(payload$mixture_path)
  scrna <- read_matrix(payload$scrna_path)
  cell_types <- as.character(read_vector(payload$cell_types_path))
  cell_states <- cell_types
  if (is_present(payload$cell_states_path)) {
    cell_states <- as.character(read_vector(payload$cell_states_path))
  }

  ref <- InstaPrism::refPrepare(
    sc_Expr = scrna,
    cell.type.labels = cell_types,
    cell.state.labels = cell_states
  )
  result <- InstaPrism::InstaPrism(
    bulk_Expr = bulk,
    refPhi_cs = ref
  )
  theta <- t(result@Post.ini.ct@theta)
  write_output(as.matrix(theta), payload$output_path)

  if (is_present(payload$extra_output_dir)) {
    dir.create(payload$extra_output_dir, showWarnings = FALSE, recursive = TRUE)
    write_df_tsv(theta, file.path(payload$extra_output_dir, "theta.tsv"))
    if (as.logical(payload$write_z)) {
      z <- InstaPrism::get_Z_array(result)
      saveRDS(z, file = file.path(payload$extra_output_dir, "Z_array.rds"))
    }
    manifest <- data.frame(
      key = if (as.logical(payload$write_z)) c("theta", "Z_array") else c("theta"),
      file = if (as.logical(payload$write_z)) c("theta.tsv", "Z_array.rds") else c("theta.tsv"),
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
} else if (tool == "quantiseq") {
  run_quantiseq(payload)
} else if (tool == "mcp_counter") {
  run_mcp_counter(payload)
} else if (tool == "psea") {
  run_psea(payload)
} else if (tool == "cellcode") {
  run_cellcode(payload)
} else if (tool == "abis") {
  run_abis(payload)
} else if (tool == "estimate") {
  run_estimate(payload)
} else if (tool == "consensus_tme") {
  run_consensus_tme(payload)
} else if (tool == "timer") {
  run_timer(payload)
} else if (tool == "cdseq") {
  run_cdseq(payload)
} else if (tool == "xcell2") {
  run_xcell2(payload)
} else if (tool == "xcell2_train") {
  run_xcell2_train(payload)
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
} else if (tool == "music2") {
  run_music2(payload)
} else if (tool == "music2_sce") {
  run_music2_sce(payload)
} else if (tool == "scdc") {
  run_scdc(payload)
} else if (tool == "mead") {
  run_mead(payload)
} else if (tool == "instaprism") {
  run_instaprism(payload)
} else if (tool == "bisque") {
  run_bisque(payload)
} else if (tool == "bayesprism") {
  run_bayesprism(payload)
} else {
  stop(paste0("Unsupported R tool: ", tool))
}
