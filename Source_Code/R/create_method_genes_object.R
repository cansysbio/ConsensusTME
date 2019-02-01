createMethodGenesList <- function(saveRdata = FALSE, path.save = '.', outputRawGenes = FALSE){
# Creates the methods genes list containing the gene sets
# for each method.
#
# Args:
#  saveRdata: If True, the method genes list will be saved as an RData file (in working directory by default)
#  path.save: When Save == TRUE a path can be specified to save the methods RData file to
#  outputRawGenes: When TRUE will output an Rdata object containing the raw gene sets from each method
#  
# Returns:
#   Returns a list containing the genes extracted from each of the methods

  require(reshape2)
  require(plyr)

  source('~/Documents/R/useful_functions.R')

  # Path for method data
  path.load = '~/Documents/ConsensusTME_Work/Consensus_Git/ConsensusTME/Source_Code/data/Method_genes/'

  #### Read all method genes in a lists ####

  # Bindea Gene Sets
  bindea.list <- as.list(read.table(sprintf('%s/Bindea_genes.txt',path.load),
                                    sep = '\t', header = T, stringsAsFactors = F))

  # Danaher Gene Sets
  danaher.list <- as.list(read.table(sprintf('%s/Danaher_genes.txt',path.load),
                                     sep = '\t', header = T, stringsAsFactors = F))
  
  # Davoli Gene Sets
  davoli.list <- as.list(read.table(sprintf('%s/Davoli_genes.txt',path.load),
                               sep = '\t', header = T, stringsAsFactors = F))

  # MCP-Counter Gene Sets
  mcp <- read.table(sprintf('%s/MCP_genes.txt',path.load),
                    sep = '\t', header = T, stringsAsFactors = F)
  mcp <- mcp[ ,colnames(mcp)!='ENTREZID']
  mcp$ind <- with(mcp, ave(seq_along(Cell.population), Cell.population, FUN=seq_along))
  mcp.list <- as.list(dcast(mcp, ind ~ Cell.population, value.var = 'HUGO.symbols',fill = '')[ ,-1])

  # xCell Gene Sets
  xcell.list <- as.list(read.table(sprintf('%s/xCell_genes.txt',path.load),
                                   sep = '\t', header = T, stringsAsFactors = F))
  names(xcell.list) <- sapply(strsplit(names(xcell.list),'_'), function(x) x [[1]])
  xcell.list <- mergeDups(xcell.list)

  # CIBERSORT Gene Sets From LM22
  # Genes with Z-Score > 1.96 are selected
  cibersort <- read.table(sprintf('%s/CIBERSORT_LM22.txt',path.load),
                           sep = '\t', header = T, row.names = 'Gene.symbol')
  cib.zscores <- apply(cibersort, 2, function(x) scale(x, center = TRUE, scale = TRUE))
  row.names(cib.zscores) <- row.names(cibersort)
  cibersort.list <- vector('list',length(colnames(cib.zscores)))
  for (i in 1:length(colnames(cib.zscores))){
    tmp.genes <- names(cib.zscores[cib.zscores[ ,i] > 1.96, i])
    cibersort.list[[i]]<- tmp.genes
  }
  names(cibersort.list) <- colnames(cib.zscores)

  # TIMER Gene Sets With Negative Purity Correlation  (All Cancers)
  timer.all <- read.table(sprintf('%s/TIMER_genes.txt',path.load),
                          sep = '\t', header = T, stringsAsFactors = F)

  # ESTIMATE Gene Sets With Negative Purity Correlation (All Cancers)
  estimate.all <- read.table(sprintf('%s/ESTIMATE_stromal_genes.txt',path.load),
                              sep = '\t', header = T, stringsAsFactors = F)

  #### Create Object with Raw Gene Sets ####
  
  raw.genes <- list()
  raw.genes$Bindea <- bindea.list
  raw.genes$CIBERSORT <- cibersort.list
  raw.genes$Davoli <- davoli.list
  raw.genes$Danaher <- danaher.list
  raw.genes$MCP.Counter <- mcp.list
  raw.genes$xCell <- xcell.list
  raw.genes$TIMER <- timer.all
  raw.genes$ESTIMATE <- estimate.all
  
  if(outputRawGenes){
    save(raw.genes, file = sprintf('%s/Raw_Gene_List.RData', path.save))
  }
  
  #### Ensure Consistant Nomenclature For Cell Types ####

  # Bindea
  bindea.prep <- list()
  bindea.prep$B_cells <- bindea.list$B_cells
  bindea.prep$Cytotoxic_cells <- bindea.list$Citotoxic_cells
  bindea.prep$Dendritic_cells <- c(bindea.list$activated_Dendritic_cells,bindea.list$immature_Dendritic_cells,
                                  bindea.list$Dendritic_cells,bindea.list$plasmacytoid_Dendritic_cells)
  bindea.prep$Eosinophils <- bindea.list$Eosinophils
  bindea.prep$Macrophages <- bindea.list$Macrophages
  bindea.prep$Mast_cells <- bindea.list$Mast_cells
  bindea.prep$NK_cells <- c(bindea.list$NK_cells,bindea.list$NK_CD56bright_cells,bindea.list$NK_CD56dim_cells)
  bindea.prep$Neutrophils <- bindea.list$Neutrophils
  bindea.prep$T_cells_CD4 <- bindea.list$T_helper_cells
  bindea.prep$T_cells_CD8 <- bindea.list$CD8_T_cells
  bindea.prep$T_cells_gamma_delta <- bindea.list$T_gamma_delta_cells
  bindea.prep$T_regulatory_cells <- bindea.list$Treg_cells

  bindea.prep <- removeBlanks(bindea.prep)
  bindea.prep <- lapply(bindea.prep, unique)

  # Danaher
  
  danaher.prep <- list()
  danaher.prep$B_cells <- danaher.list$B.cells
  danaher.prep$Cytotoxic_cells <- danaher.list$Cytotoxic.cells
  danaher.prep$Dendritic_cells <- danaher.list$DC
  danaher.prep$Macrophages <- danaher.list$Macrophages
  danaher.prep$Mast_cells <- danaher.list$Mast.cells
  danaher.prep$Neutrophils <- danaher.list$Neutrophils
  danaher.prep$NK_cells <- c(danaher.list$NK.CD56dim.cells, danaher.list$NK.cells)
  danaher.prep$T_cells_CD4 <- danaher.list$Th1.cells
  danaher.prep$T_cells_CD8 <- c(danaher.list$CD8.T.cells, danaher.list$Exhausted.CD8)
  danaher.prep$T_regulatory_cells <- danaher.list$Treg
  
  danaher.prep <- removeBlanks(danaher.prep)
  danaher.prep <- lapply(danaher.prep, unique)
  
  # Davoli

  davoli.prep <- list()
  davoli.prep$B_cells <- davoli.list$B_cells
  davoli.prep$Dendritic_cells <- davoli.list$Dendritics
  davoli.prep$Macrophages <- davoli.list$Macrophages
  davoli.prep$Macrophages_M1 <- davoli.list$Macrophages_M1
  davoli.prep$Macrophages_M2 <- davoli.list$Macrophages_M2
  davoli.prep$NK_cells <- c(davoli.list$NK_cells, davoli.list$CD8_effector_NK_cells)
  davoli.prep$T_cells_CD4 <- davoli.list$CD4_mature
  davoli.prep$T_cells_CD8 <- davoli.list$CD8_effector
  davoli.prep$T_regulatory_cells <- davoli.list$T_regs

  davoli.prep <- removeBlanks(davoli.prep)
  davoli.prep <- lapply(davoli.prep, unique)

  # MCP-Counter

  mcp.prep <- list()
  mcp.prep$B_cells <- mcp.list$`B lineage`
  mcp.prep$Cytotoxic_cells <- mcp.list$`Cytotoxic lymphocytes`
  mcp.prep$Dendritic_cells <- mcp.list$`Myeloid dendritic cells`
  mcp.prep$Endothelial <- mcp.list$`Endothelial cells`
  mcp.prep$Fibroblasts <- mcp.list$Fibroblasts
  mcp.prep$NK_cells <- mcp.list$`NK cells`
  mcp.prep$Neutrophils <- mcp.list$Neutrophils
  mcp.prep$T_cells_CD8 <- mcp.list$`CD8 T cells`

  mcp.prep <- removeBlanks(mcp.prep)
  mcp.prep <- lapply(mcp.prep, unique)

  # xCell

  xcell.prep <- list()
  xcell.prep$B_cells <- c(xcell.list$B.cells, xcell.list$Class.switched.memory.B.cells,
                          xcell.list$Memory.B.cells, xcell.list$naive.B.cells, xcell.list$pro.B.cells)
  xcell.prep$Dendritic_cells <- c(xcell.list$DC, xcell.list$aDC, xcell.list$cDC,
                                  xcell.list$iDC, xcell.list$pDC)
  xcell.prep$Endothelial <- c(xcell.list$Endothelial.cells, xcell.list$ly.Endothelial.cells, xcell.list$mv.Endothelial.cells)
  xcell.prep$Eosinophils <- xcell.list$Eosinophils
  xcell.prep$Fibroblasts <- xcell.list$Fibroblasts
  xcell.prep$Macrophages <- c(xcell.list$Macrophages, xcell.list$Macrophages.M1, xcell.list$Macrophages.M2)
  xcell.prep$Macrophages_M1 <- xcell.list$Macrophages.M1
  xcell.prep$Macrophages_M2 <- xcell.list$Macrophages.M2
  xcell.prep$Mast_cells <- xcell.list$Mast.cells
  xcell.prep$Monocytes <- xcell.list$Monocytes
  xcell.prep$NK_cells <- xcell.list$NK.cells
  xcell.prep$Neutrophils <- xcell.list$Neutrophils
  xcell.prep$Plasma_cells <- xcell.list$Plasma.cells
  xcell.prep$T_cells_CD4 <- c(xcell.list$CD4..T.cells, xcell.list$CD4..Tcm, xcell.list$CD4..Tem,
                              xcell.list$CD4..memory.T.cells, xcell.list$CD4..naive.T.cells)
  xcell.prep$T_cells_CD8 <- c(xcell.list$CD8..T.cells, xcell.list$CD8..Tcm, xcell.list$CD8..Tem,
                              xcell.list$CD8..naive.T.cells)
  xcell.prep$T_cells_gamma_delta <- xcell.list$Tgd.cells
  xcell.prep$T_regulatory_cells <- xcell.list$Tregs

  xcell.prep <- removeBlanks(xcell.prep)
  xcell.prep <- lapply(xcell.prep, unique)

  # CIBERSORT

  cibersort.prep <- list()
  cibersort.prep$B_cells <- c(cibersort.list$B.cells.memory, cibersort.list$B.cells.naive)
  cibersort.prep$Dendritic_cells <- c(cibersort.list$Dendritic.cells.activated, cibersort.list$Dendritic.cells.resting)
  cibersort.prep$Eosinophils <- cibersort.list$Eosinophils
  cibersort.prep$Macrophages <- cibersort.list$Macrophages.M0
  cibersort.prep$Macrophages_M1 <- cibersort.list$Macrophages.M1
  cibersort.prep$Macrophages_M2 <- cibersort.list$Macrophages.M2
  cibersort.prep$Mast_cells <- c(cibersort.list$Mast.cells.activated, cibersort.list$Mast.cells.resting)
  cibersort.prep$Monocytes <- cibersort.list$Monocytes
  cibersort.prep$NK_cells <- c(cibersort.list$NK.cells.activated, cibersort.list$NK.cells.resting)
  cibersort.prep$Neutrophils <- cibersort.list$Neutrophils
  cibersort.prep$Plasma_cells <- cibersort.list$Plasma.cells
  cibersort.prep$T_cells_CD4 <- c(cibersort.list$T.cells.CD4.memory.activated, cibersort.list$T.cells.CD4.memory.resting,
                                  cibersort.list$T.cells.CD4.naive)
  cibersort.prep$T_cells_CD8 <- cibersort.list$T.cells.CD8
  cibersort.prep$T_cells_gamma_delta <- cibersort.list$T.cells.gamma.delta
  cibersort.prep$T_regulatory_cells <- cibersort.list$T.cells.regulatory..Tregs.

  cibersort.prep <- removeBlanks(cibersort.prep)
  cibersort.prep <- lapply(cibersort.prep, unique)

  #### Create Method Data Object ####

  method.genes <- list()
  method.genes$Bindea <- bindea.prep
  method.genes$CIBERSORT <- cibersort.prep
  method.genes$Davoli <- davoli.prep
  method.genes$Danaher <- danaher.prep
  method.genes$MCP.Counter <- mcp.prep
  method.genes$xCell <- xcell.prep
  method.genes$TIMER <- timer.all
  method.genes$ESTIMATE <- estimate.all

  #### Save Method Data Object ####
  if(saveRdata){
    save(method.genes, file = sprintf('%s/Method_data.RData', path.save))
  }
  return(method.genes)
}
