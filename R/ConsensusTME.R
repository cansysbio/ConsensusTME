#' ConsensusTME: A package for estimation of cell type abundance.
#'
#' @docType package
#' @name ConsensusTME
NULL

#' Raw Method Signatures and Filter Data
#'
#' Compilation of gene sets, signature matracies and filters for generation
#' of ConsensusTME gene sets
#' @docType data
#'
#' @usage data(methodSignatures)
#'
#' @format list:
#' \describe{
#'   \item{Bindea}{Gene sets for 24 immune cells & 2 biological processes}
#'   \item{Danaher}{Gene sets for 14 immune cells}
#'   \item{Davoli}{Gene sets for 10 immune cells}
#'   \item{MCP.Counter}{Gene sets for 10 immune cells & 2 stromal cells}
#'   \item{xCell}{Gene sets for 34 immune cells, 9 other cells of haematopoietic lineage &
#'   21 cells of non haematopoietic lineage}
#'   \item{CIBERSORT}{LM22 signature matrix for 22 immune cell types}
#'   \item{ImmuneGenes}{Dataframe with immune genes with negative correlation with tumour
#'   purity for each cancer type}
#'   \item{ESTIMATE}{Dataframe with ESTIMATES stromal genes that have a negative correlation
#'   with tumour purity for each cancer type}
#' }
#' @export
"methodSignatures"

#' Consensus TME Cancer Types
#'
#' TCGA cancer types for which ConsensusTME can generate gene sets for:
#'
#' ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC LUAD LUSC MESO OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM.
#' @docType data
#'
#' @usage data(cancerAll)
#'
#' @description character vector containing TCGA abbreviations
#'
#' @export
"cancerAll"

#' Consensus TME Gene Sets
#'
#' Preprocessed ConsensusTME gene sets for each of the TCGA cancer types
#'
#' @docType data
#'
#' @usage data(consensusGeneSets)
#'
#' @description list of genes for each cell type for each cancer type
#'
#' @export
"consensusGeneSets"

#' Run ConsensusTME Cell Type Estimation
#' \code{consensusTMEAnalysis} takes bulk tumour gene expression data and returns
#' cell type specific enrichment scores for each sample
#'
#' @param bulkExp bulk tumour gene expression matrix. HUGO gene symbols as row names & sample IDs
#' as column names.
#' @param cancerType string passed to indicate which TCGA cancer type samples are most similar to. \bold{N.B} samples of different cancer types should be run seperately.
#' Available cancer types: \code{"ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP","LGG", "LIHC", "LUAD", "LUSC", "MESO",
#' "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "Unfiltered"}.
#' @param statMethod statistical framework to be used in generating gene set enrichment scores.
#' These mirror the parameter options of \code{GSVA::gsva()} with the exception of \code{singScore}.
#' which leverages \code{singscore::multiScore()}. Default: \code{ssgsea}
#' @param singScoreDisp logical, when using singscore method should the dispersion for each
#' gene set be returned with the enrichment scores. Default: \code{FALSE}
#' @param immuneScore logical, when \code{TRUE} (default) an Immune Score is produced representing overall
#' level of immune infiltration for each sample.
#' @param excludeCells a character vector that defines cell types to not generate enrichment scores for.
#' @param parallel.sz parameter passed to \code{GSVA::gsva()} - Number of processors to use when doing the calculations in parallel.
#' This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz = 0)
#' then it will use all available core processors unless this is set with a smaller number.
#' If snow is loaded then this must be set to a positive integer number that specifies the number of processors to employ in the parallel calculation.
#'
#' @return returns estimation of cell type abundance for each sample in the bulk tumour gene expression matrix
#' @export

consensusTMEAnalysis <- function(bulkExp, cancerType = NULL, statMethod = c("ssgsea", "gsva", "plage","zscore", "singScore"),
                         singScoreDisp = FALSE, immuneScore = TRUE, excludeCells = NULL,
                         parallel.sz = 0) {
  
  if(is.null(cancerType)) {
    stop(paste0("argument \"cancerType\" is missing and should be one of: ", paste0(cancerAll, collapse = ", "), ", Unfiltered"))
  }
  cancerType <- match.arg(cancerType, c(cancerAll, "Unfiltered"))
  statMethod <- match.arg(statMethod)
  sampleSize <- dim(bulkExp)[2]

  # Create Consensus Gene Sets

  matchedSigs <- matchGeneSigs(methodSignatures = methodSignatures)

  purityGenes <- methodSignatures$ImmuneGenes
  estimateGenes <- methodSignatures$ESTIMATE

  consensusGenesAll <- buildConsensusGenes(matchedSigs = matchedSigs, consCancerType = cancerType,
                                           immuneFilter = purityGenes, stromalFilter = estimateGenes)

  consensusGenes <- consensusGenesAll[[cancerType]]

  # Exclude Designated Cell Types

  consensusGenes <- consensusGenes[!names(consensusGenes) %in% excludeCells]

  # Give Warning For Cell Types Unable To Be Estimated

  sigLengths <- sapply(consensusGenes, length)

  if (any(sigLengths == 0)) {
    missingCells <- names(sigLengths[sigLengths == 0])
    consensusGenes <- consensusGenes[!sigLengths == 0]
    warning(paste0("the following signatures did not have any genes that meet the purity filter: ",
                   paste(missingCells, collapse = ", ")))

  }

  cat(paste0("Producing ConsensusTME Estimates Using The Following Parameters:",
             "\n Statistical Framework: \"", statMethod, "\"",
             "\n Gene Sets For Cancer Type: \"", cancerType, "\"",
             "\n Sample Size: ", sampleSize, "\n"))

  # Run Estimation

  consensusScores <- geneSetEnrichment(bulkExp = bulkExp, signatures = consensusGenes, statMethod = statMethod, singScoreDisp = singScoreDisp)

  return(consensusScores)

}

#' Carry out Gene Set Enrichment Analysis
#'
#' \code{geneSetEnrichment} Combines exisiting signatures and allows additional of a new signature.
#'
#' @param bulkExp bulk tumour gene expression matrix. HUGO gene symbols as row names & sample IDs
#' as column names.
#' @param signatures a list with each element containing genes to represent a cell type. The cell types
#' should be the names of each element of the list.
#' @param statMethod statistical framework to be used in generating gene set enrichment scores.
#' These mirror the parameter options of \code{GSVA::gsva()} with the exception of \code{singScore}.
#' which leverages \code{singscore::multiScore()}. Default: \code{ssgsea}
#' @param singScoreDisp Logical, when using singscore method should the dispersion for each
#' gene set be returned with the enrichment scores. Default: \code{FALSE}
#' @param parallel.sz parameter for \code{GSVA::gsva()} - Number of processors to use when doing the calculations in parallel.
#' This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz = 0)
#' then it will use all available core processors unless this is set with a smaller number.
#' If snow is loaded then this must be set to a positive integer number that specifies the number of processors to employ in the parallel calculation.
#'
#' @return returns a list with curated signatures ready to be combined
#' @export


geneSetEnrichment <- function(bulkExp, signatures, statMethod = c("ssgsea", "gsva", "plage", "zscore", "singScore"), singScoreDisp = FALSE,
                              parallel.sz = 0) {

  statMethod <- match.arg(statMethod)

  # Check for numerical matrix and ensure there's enough genes to run gene set enrichment

  if (!(is.matrix(bulkExp) & is.numeric(bulkExp))) {
    stop("Error: bulkExp must be a numerical matrix")
  } else if (dim(bulkExp)[1] < 1000) {
    warning("Low number of genes in bulk expression matrix, enrichment analysis may not be reliable")
  }

  # Ensure Overlap With Bulk Matrix Genes

  bulkGenes <- row.names(bulkExp)
  intLengths <- sapply(signatures, function(x) length(intersect(bulkGenes, x)))

  if (all(intLengths == 0)) {
    stop("Error: No genes in the gene sets could be matched to the identifiers in the expression data")
  } else if (any(intLengths == 0)) {
    depletedCells <- names(intLengths[intLengths == 0])
    signatures <- signatures[!names(signatures) %in% depletedCells]
    warning(paste0("The following signatures will not be ran as no genes could be matched to the identifiers in the expression data: ",
                   paste(depletedCells, collapse = ", ")))
  }

  # Estimation Methods

  if (statMethod == "singScore") {
    rankedExp <- singscore::rankGenes(bulkExp)
    signatures <- signatures[sapply(signatures, function(x) length(x) > 0)]
    sigsGeneSet <- lapply(names(signatures), function(cellType){
      GSEABase::GeneSet(signatures[[cellType]], setName = cellType)
    })
    sigsGeneCol <- GSEABase::GeneSetCollection(sigsGeneSet)
    singOut <- singscore::multiScore(rankData = rankedExp, upSetColc = sigsGeneCol)
    if (!singScoreDisp) {
      singOut$Dispersions <- NULL
    }
    consensusEstimates <- singOut
  } else if (statMethod %in% c("ssgsea", "gsva", "plage", "zscore")) {
    if(packageVersion("GSVA") >= "1.50.0"){
      # Legacy GSVA function was depreciated in version 1.50 and entirely removed in version 1.52. (minSize must be 1 or greater)
      if (statMethod=="ssgsea") {
        gsvapar <- GSVA::ssgseaParam(exprData = bulkExp, geneSets = signatures, minSize = 1, maxSize = dim(bulkExp)[1], normalize = TRUE)
      } else if (statMethod=="gsva") {
        gsvapar <- GSVA::gsvaParam(exprData = bulkExp, geneSets = signatures, minSize = 1, maxSize = dim(bulkExp)[1])
      } else if (statMethod=="plage") {
        gsvapar <- GSVA::plageParam(exprData = bulkExp, geneSets = signatures, minSize = 1, maxSize = dim(bulkExp)[1])
      } else {
        gsvapar <- GSVA::zscoreParam(exprData = bulkExp, geneSets = signatures, minSize = 1, maxSize = dim(bulkExp)[1])
      }
      consensusEstimates <- GSVA::gsva(gsvapar)
    } else {
      consensusEstimates <- GSVA::gsva(expr = bulkExp, gset.idx.list = signatures, method = statMethod, min.sz = 0, max.sz = dim(bulkExp)[1],
                                       parallel.sz = parallel.sz, ssgsea.norm = TRUE)
    }
  }
}

#' Build Consensus Gene Supersets
#' \code{buildConsensusGenes} Combines pre-processed signatures and filters to create curated gene sets.
#'
#' @param matchedSigs list of pre-processed signatures with consistant cell type nomenclature
#' @param consCancerType TCGA cancer types to produce gene sets for if no cancer type is specified Default: All TCGA cancer types.
#' @param immuneFilter list of genes to filter immune gene sets by. If \code{NULL} - \code{methodSignatures$ImmuneGenes}
#' is used (list of immune genes with negative correlation with tumour purity for each tcga cancer type).
#' \bold{N.B.} Must have two columns: "Gene_Symbol" & "Cancer".
#' @param stromalFilter list of genes to filter stromal gene sets by. If \code{NULL} - \code{methodSignatures$ESTIMATE}
#' is used (list of ESTIMATE's stromal genes with negative correlation with tumour purity for each tcga cancer type).
#' \bold{N.B.} Must have two columns: "Gene_Symbol" & "Cancer".
#' @param immuneScore logical, when \code{TRUE} (default) an Immune Score is produced representing overall
#' level of immune infiltration for each sample.
#' @param unfilteredGeneSet logical, include geneset without cancer specific filtering. If left empty & more than one cancer type has been
#' specified defaults to \code{TRUE} otherwise defaults to \code{FALSE}.
#'
#' @return returns consensusGeneSets
#' @export

buildConsensusGenes <- function(matchedSigs, consCancerType = c(cancerAll, "Unfiltered"), immuneFilter = NULL, stromalFilter = NULL, immuneScore = TRUE,
                                unfilteredGeneSet = NULL) {

  consCancerType <- match.arg(consCancerType, c(cancerAll, "Unfiltered"), several.ok = TRUE)

  if (is.null(immuneFilter)) {
    immuneFilter <- methodSignatures$ImmuneGenes
  }
  if (is.null(stromalFilter)) {
    stromalFilter <- methodSignatures$ESTIMATE
  }

  # Merge Genes For Each Cell Type Across All Methods

  geneSigs <- matchedSigs[!(names(matchedSigs) %in% c("ImmuneGenes", "ESTIMATE"))]

  cellTypes <- unique(unlist(lapply(geneSigs, names)))
  consensusAll <- lapply(cellTypes, function(cell){
    allGenes <- lapply(geneSigs, function (x) x[names(x) == cell])
    return(unique(unlist(allGenes)))
  })
  names(consensusAll) <- cellTypes

  # Apply Gene Filters, Creating Cancer Specific Gene Sets

  stromalCells <- c("Endothelial", "Fibroblasts")
  immuneCells <- cellTypes[!cellTypes %in% stromalCells]

  consensusGenes <- lapply(consCancerType, function(canc) {
    if (canc == "Unfiltered") {
      if (immuneScore) {
        consensusAll$Immune_Score <- sort(unique(unlist(consensusAll)))
      }
      return(consensusAll)
    }
    immuneGenes <- immuneFilter[immuneFilter$Cancer == canc, 1]
    stromalGenes <- stromalFilter[stromalFilter$Cancer == canc, 1]
    filteredSigs <- lapply(cellTypes, function(cell){
      if (cell %in% immuneCells) {
        return(sort(intersect(consensusAll[[cell]], immuneGenes)))
      } else if (cell %in% stromalCells) {
        return(sort(intersect(consensusAll[[cell]], stromalGenes)))
      } else {
        stop("Error: cell type not in stromal or immune cell types")
      }
    })
    names(filteredSigs) <- cellTypes

    if (immuneScore) {
      immSigs <- filteredSigs[!names(filteredSigs) %in% stromalCells]
      filteredSigs$Immune_Score <- sort(unique(unlist(immSigs)))
    }

    return(filteredSigs)
  })

  names(consensusGenes) <- consCancerType

  # Add Unfiltered Set of Genes If Required

  if (is.null(unfilteredGeneSet) & length(consCancerType) > 1 & !"Unfiltered" %in% names(consensusGenes)) {
    consensusGenes$Unfiltered <- consensusAll
  }


  return(consensusGenes)

}

#' Combine Existing Signatures Ensuring Consistant Nomenclature
#'
#' \code{matchGeneSigs} Combines existing signatures.
#'
#' @param methodSignatures raw gene signatures list generated from methodSignatures.R
#'
#' @return returns a list with curated signatures ready to be combined
#' @export


matchGeneSigs <- function(methodSignatures) {

  matchedSigs <- list()

  # Bindea

  bindeaMatched <- list()

  bindeaRaw <- methodSignatures$Bindea

  bindeaMatched$B_cells <- bindeaRaw$B_cells
  bindeaMatched$Cytotoxic_cells <- bindeaRaw$Citotoxic_cells
  bindeaMatched$Dendritic_cells <- c(bindeaRaw$activated_Dendritic_cells,
                                     bindeaRaw$immature_Dendritic_cells,
                                     bindeaRaw$Dendritic_cells,
                                     bindeaRaw$plasmacytoid_Dendritic_cells)
  bindeaMatched$Eosinophils <- bindeaRaw$Eosinophils
  bindeaMatched$Macrophages <- bindeaRaw$Macrophages
  bindeaMatched$Mast_cells <- bindeaRaw$Mast_cells
  bindeaMatched$NK_cells <- c(bindeaRaw$NK_cells,
                              bindeaRaw$NK_CD56bright_cells,
                              bindeaRaw$NK_CD56dim_cells)
  bindeaMatched$Neutrophils <- bindeaRaw$Neutrophils
  bindeaMatched$T_cells_CD4 <- bindeaRaw$T_helper_cells
  bindeaMatched$T_cells_CD8 <- bindeaRaw$CD8_T_cells
  bindeaMatched$T_cells_gamma_delta <- bindeaRaw$T_gamma_delta_cells
  bindeaMatched$T_regulatory_cells <- bindeaRaw$Treg_cells

  matchedSigs$Bindea <- lapply(bindeaMatched, unique)

  # Danaher

  danaherMatched <- list()

  danaherRaw <- methodSignatures$Danaher

  danaherMatched$B_cells <- danaherRaw$B.cells
  danaherMatched$Cytotoxic_cells <- danaherRaw$Cytotoxic.cells
  danaherMatched$Dendritic_cells <- danaherRaw$DC
  danaherMatched$Macrophages <- danaherRaw$Macrophages
  danaherMatched$Mast_cells <- danaherRaw$Mast.cells
  danaherMatched$Neutrophils <- danaherRaw$Neutrophils
  danaherMatched$NK_cells <- c(danaherRaw$NK.CD56dim.cells,
                               danaherRaw$NK.cells)
  danaherMatched$T_cells_CD4 <- danaherRaw$Th1.cells
  danaherMatched$T_cells_CD8 <- c(danaherRaw$CD8.T.cells,
                                  danaherRaw$Exhausted.CD8)
  danaherMatched$T_regulatory_cells <- danaherRaw$Treg

  matchedSigs$Danaher <- lapply(danaherMatched, unique)

  # Davoli

  davoliMatched <- list()

  davoliRaw <- methodSignatures$Davoli

  davoliMatched$B_cells <- davoliRaw$B_cells
  davoliMatched$Dendritic_cells <- davoliRaw$Dendritics
  davoliMatched$Macrophages <- davoliRaw$Macrophages
  davoliMatched$Macrophages_M1 <- davoliRaw$Macrophages_M1
  davoliMatched$Macrophages_M2 <- davoliRaw$Macrophages_M2
  davoliMatched$NK_cells <- c(davoliRaw$NK_cells,
                              davoliRaw$CD8_effector_NK_cells)
  davoliMatched$T_cells_CD4 <- davoliRaw$CD4_mature
  davoliMatched$T_cells_CD8 <- davoliRaw$CD8_effector
  davoliMatched$T_regulatory_cells <- davoliRaw$T_regs

  matchedSigs$Davoli <- lapply(davoliMatched, unique)

  # MCP-Counter

  mcpMatched <- list()

  mcpRaw <- methodSignatures$MCP.Counter

  mcpMatched$B_cells <- mcpRaw$`B lineage`
  mcpMatched$Cytotoxic_cells <- mcpRaw$`Cytotoxic lymphocytes`
  mcpMatched$Dendritic_cells <- mcpRaw$`Myeloid dendritic cells`
  mcpMatched$Endothelial <- mcpRaw$`Endothelial cells`
  mcpMatched$Fibroblasts <- mcpRaw$Fibroblasts
  mcpMatched$NK_cells <- mcpRaw$`NK cells`
  mcpMatched$Neutrophils <- mcpRaw$Neutrophils
  mcpMatched$T_cells_CD8 <- mcpRaw$`CD8 T cells`

  matchedSigs$MCP.Counter <- lapply(mcpMatched, unique)

  # xCell

  xCellMatched <- list()

  xCellRaw <- methodSignatures$xCell

  xCellMatched$B_cells <- c(xCellRaw$B.cells,
                            xCellRaw$Class.switched.memory.B.cells,
                            xCellRaw$Memory.B.cells,
                            xCellRaw$naive.B.cells,
                            xCellRaw$pro.B.cells)
  xCellMatched$Dendritic_cells <- c(xCellRaw$DC,
                                    xCellRaw$aDC,
                                    xCellRaw$iDC,
                                    xCellRaw$pDC)
  xCellMatched$Endothelial <- c(xCellRaw$Endothelial.cells,
                                xCellRaw$ly.Endothelial.cells,
                                xCellRaw$mv.Endothelial.cells)
  xCellMatched$Eosinophils <- xCellRaw$Eosinophils
  xCellMatched$Fibroblasts <- xCellRaw$Fibroblasts
  xCellMatched$Macrophages <- c(xCellRaw$Macrophages,
                                xCellRaw$Macrophages.M1,
                                xCellRaw$Macrophages.M2)
  xCellMatched$Macrophages_M1 <- xCellRaw$Macrophages.M1
  xCellMatched$Macrophages_M2 <- xCellRaw$Macrophages.M2
  xCellMatched$Mast_cells <- xCellRaw$Mast.cells
  xCellMatched$Monocytes <- xCellRaw$Monocytes
  xCellMatched$NK_cells <- xCellRaw$NK.cells
  xCellMatched$Neutrophils <- xCellRaw$Neutrophils
  xCellMatched$Plasma_cells <- xCellRaw$Plasma.cells
  xCellMatched$T_cells_CD4 <- c(xCellRaw$CD4..T.cells,
                                xCellRaw$CD4..Tcm,
                                xCellRaw$CD4..Tem,
                                xCellRaw$CD4..memory.T.cells,
                                xCellRaw$CD4..naive.T.cells)
  xCellMatched$T_cells_CD8 <- c(xCellRaw$CD8..T.cells,
                                xCellRaw$CD8..Tcm,
                                xCellRaw$CD8..Tem,
                                xCellRaw$CD8..naive.T.cells)
  xCellMatched$T_cells_gamma_delta <- xCellRaw$Tgd.cells
  xCellMatched$T_regulatory_cells <- xCellRaw$Tregs

  matchedSigs$xCell <- lapply(xCellMatched, unique)

  # CIBERSORT

  # Take features from LM22 Signature Matrix
  # Genes With Z Score >1.96 are selected for each cell type

  cibersortMat <- methodSignatures$CIBERSORT

  cibersortRaw <- apply(cibersortMat, 2, function(cell){
    cellNorm <- scale(cell, center = TRUE, scale = TRUE)
    row.names(cellNorm) <- row.names(cibersortMat)
    names(cellNorm[cellNorm > 1.96, ])
  })

  names(cibersortRaw) <- colnames(cibersortMat)

  # Match Names

  cibersortMatched <- list()
  cibersortMatched$B_cells <- c(cibersortRaw$B.cells.memory,
                                cibersortRaw$B.cells.naive)
  cibersortMatched$Dendritic_cells <- c(cibersortRaw$Dendritic.cells.activated,
                                        cibersortRaw$Dendritic.cells.resting)
  cibersortMatched$Eosinophils <- cibersortRaw$Eosinophils
  cibersortMatched$Macrophages <- cibersortRaw$Macrophages.M0
  cibersortMatched$Macrophages_M1 <- cibersortRaw$Macrophages.M1
  cibersortMatched$Macrophages_M2 <- cibersortRaw$Macrophages.M2
  cibersortMatched$Mast_cells <- c(cibersortRaw$Mast.cells.activated,
                                   cibersortRaw$Mast.cells.resting)
  cibersortMatched$Monocytes <- cibersortRaw$Monocytes
  cibersortMatched$NK_cells <- c(cibersortRaw$NK.cells.activated,
                                 cibersortRaw$NK.cells.resting)
  cibersortMatched$Neutrophils <- cibersortRaw$Neutrophils
  cibersortMatched$Plasma_cells <- cibersortRaw$Plasma.cells
  cibersortMatched$T_cells_CD4 <- c(cibersortRaw$T.cells.CD4.memory.activated,
                                    cibersortRaw$T.cells.CD4.memory.resting,
                                    cibersortRaw$T.cells.CD4.naive)
  cibersortMatched$T_cells_CD8 <- cibersortRaw$T.cells.CD8
  cibersortMatched$T_cells_gamma_delta <- cibersortRaw$T.cells.gamma.delta
  cibersortMatched$T_regulatory_cells <- cibersortRaw$T.cells.regulatory..Tregs.

  matchedSigs$CIBERSORT <- lapply(cibersortMatched, unique)

  return(matchedSigs)
}

.onLoad <- function(libname, pkgname) {
  utils::data("methodSignatures", package = pkgname, envir = parent.env(environment()))
  utils::data("cancerAll", package = pkgname, envir = parent.env(environment()))
  utils::data("consensusGeneSets", package = pkgname, envir = parent.env(environment()))
}
