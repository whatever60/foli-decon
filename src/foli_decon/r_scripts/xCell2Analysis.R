xCell2Analysis <- function(mix,
                           xcell2object,
                           minSharedGenes = 0.9,
                           rawScores = FALSE,
                           spillover = TRUE,
                           spilloverAlpha = 0.5,
                           BPPARAM = BiocParallel::SerialParam()) {
  
  message("Starting xCell2 Analysis...")
  
  scoreMix <- function(cellType, mixRanked, signaturesCellType) {

    scores <- vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(mixRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(mixRanked)))
    scores <- matrix(scores, nrow = ncol(mixRanked))
    rownames(scores) <- colnames(mixRanked)
    return(scores)
  }
  calcEnrichment <- function(cellType) {
    
    signaturesCellType <- getSignatures(xcell2object)[startsWith(names(
      getSignatures(xcell2object)
    ), paste0(cellType, "#"))]
    
    # Remove signature with missing genes
    sigs2remove <- vapply(signaturesCellType, function(sig){
      sum(sig %in% rownames(mixRanked)) == 0
    }, FUN.VALUE = logical(1))
    
    if (all(sigs2remove) | sum(!sigs2remove) < 3) {
      warning(
        "Cannot calculate enrichment scores for '",
        cellType, "' because all genes in those signatures are missing in your mixture."
      )
      zeros_matrix <- matrix(rep(0, ncol(mixRanked)*3), nrow = ncol(mixRanked), ncol = 3)
      rownames(zeros_matrix) <- colnames(mixRanked)
      return(zeros_matrix)
    }
    
    if (any(sigs2remove)) {
      warning(
        "Removing ", sum(sigs2remove), "/", length(signaturesCellType), " signatures of '",
        cellType, "' because all genes in those signatures are missing in your mixture."
      )
    }
    
    signaturesCellType <- signaturesCellType[!sigs2remove]
    
    # Calculate enrichment scores with SingScore
    warning_list <- list()
    
    scores <- withCallingHandlers(
      {
        scoreMix(cellType, mixRanked, signaturesCellType)
      },
      warning = function(w) {
        warning_list <- c(warning_list, conditionMessage(w))  # Capture warnings
        invokeRestart("muffleWarning")  # Prevent the warning from being printed repeatedly
      }
    )
    
    if (length(warning_list) > 0) {
      warning("There were ",  length(warning_list), "/", length(signaturesCellType),
              " signatures with missing gene(s) for '", cellType,
              "' during enrichment scores calculation.\n",
              "Missing genes are:\n",
              unique(unlist(strsplit(gsub(".*genes missing: ", "", unlist(warning_list)), ","))))
    }
    
    return(scores)
  }
  
  # Wrapper function to include progress bar
  calcEnrichmentWrapper <- function(cellType) {
    pb$tick()
    calcEnrichment(cellType)
  }
  
  # Check for missing values
  if (any(is.na(mix))) {
    stop("The mixture contains NA or NaN values. Please clean your data.")
  }

  # Check reference/mixture genes intersection
  shared_genes <- intersect(rownames(mix), getGenesUsed(xcell2object))
  genesIntersectFrac <- round(length(shared_genes) / length(getGenesUsed(xcell2object)), 2)
  if (genesIntersectFrac < minSharedGenes) {
    stop(
      "This xCell2 reference shares ",
      genesIntersectFrac, " genes with the mixtures and minSharedGenes = ",
      minSharedGenes, ".",
      "\n", "Consider training a new xCell2 reference or adjusting minSharedGenes."
    )
  }

  if (genesIntersectFrac < 0.85) {
    warning(
      "This xCell2 reference shares only ",
      genesIntersectFrac, " genes with the mixtures.",
      "\n", "Consider using a reference that share at least 85% of the genes with the mixture for optimal results."
    )    
  }
  
  # Rank mix gene expression matrix
  mixRanked <- singscore::rankGenes(mix[shared_genes, , drop = FALSE], tiesMethod="average")
  
  # Score and predict
  sigsCellTypes <- unique(unlist(lapply(
    names(getSignatures(xcell2object)),
    function(x) {
      strsplit(x, "#")[[1]][1]
    }
  )))
  
  # Get raw enrichment scores
  message("Calculating enrichment scores for all cell types...")
  
  if (BiocParallel::bpworkers(BPPARAM) > 1) {
    resRaw <- BiocParallel::bplapply(sigsCellTypes, calcEnrichment, BPPARAM = BPPARAM)
  }else{
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent (:current/:total) :elapsed",
      total = length(sigsCellTypes), clear = FALSE, width = 60
    )
    resRaw <- BiocParallel::bplapply(sigsCellTypes, calcEnrichmentWrapper, BPPARAM = BPPARAM)
  }

  names(resRaw) <- sigsCellTypes
  
  res <- do.call(rbind, lapply(resRaw, rowMeans))
  rownames(res) <- names(resRaw)
  colnames(res) <- rownames(resRaw[[1]])
  
  # Check for negative enrichment scores
  neg_enrichment <- names(which(apply(res, 2, function(x){any(x<0)})))
  if (length(neg_enrichment) > 0) {
    error_msg <- paste0("There is a problem with the following samples: ",
                        paste(neg_enrichment, collapse = ", "),
                        ". Please check your input.")
    stop(error_msg)
  }
  
  if (rawScores) {
    message("Returning raw enrichment scores without linear transformation or correction.")
    return(res)
  } else {
    # Linear transformation
    res <- t(vapply(rownames(res), function(cellType) {
      cellTypeRes <- res[cellType, ]
      
      # Linear transformation
      a <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$a
      b <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$b
      m <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$m
      
      cellTypeRes <- (cellTypeRes^(1 / b)) / a
      cellTypeRes <- cellTypeRes * m
      
      # Shift values
      cellTypeRes <- cellTypeRes - min(cellTypeRes)
      cellTypeRes <- round(cellTypeRes, 5)
      
      return(cellTypeRes)
    }, FUN.VALUE = double(ncol(res))))
  }
  
  if (spillover) {
    message("Performing spillover correction...")
    # Spillover correction
    spillMat <- getSpillMat(xcell2object) * spilloverAlpha
    diag(spillMat) <- 1
    
    rows <- intersect(rownames(res), rownames(spillMat))
    
    scoresCorrected <- apply(res[rows, ], 2, function(x) {
      pracma::lsqlincon(spillMat[rows, rows],
                        x,
                        lb = 0
      )
    })
    scoresCorrected[scoresCorrected < 0] <- 0
    rownames(scoresCorrected) <- rows
    
    message("xCell2 Analysis completed successfully.")
    return(scoresCorrected)
  } else {
    message("Skipping spillover correction...")
    message("xCell2 Analysis completed successfully.")
    return(res)
  }
}
