#### Create Method Signatures Objects ####

library(reshape2)

source("./data-raw/useful_functions.R")

# Initialise list

methodSignatures <- list()

# Bindea Gene Sets
methodSignatures$Bindea <- as.list(read.table("./data-raw/Bindea_genes.txt",
                                      sep = "\t",
                                      header = T,
                                      stringsAsFactors = F))

# Danaher Gene Sets
methodSignatures$Danaher <- as.list(read.table("./data-raw/Danaher_genes.txt",
                                       sep = "\t",
                                       header = T,
                                       stringsAsFactors = F))

# Davoli Gene Sets
methodSignatures$Davoli <- as.list(read.table("./data-raw/Davoli_genes.txt",
                                      sep = "\t",
                                      header = T,
                                      stringsAsFactors = F))

# MCP-Counter Gene Sets
mcp <- read.table("./data-raw/MCP_genes.txt",
                  sep = "\t",
                  header = T,
                  stringsAsFactors = F)

mcp <- mcp[ ,colnames(mcp)!="ENTREZID"]

mcp$ind <- with(mcp, ave(seq_along(Cell.population),
                         Cell.population,
                         FUN = seq_along))

methodSignatures$MCP.Counter <- as.list(dcast(mcp,
                                      ind ~ Cell.population,
                                      value.var = "HUGO.symbols",
                                      fill = "")[ ,-1])

# xCell Gene Sets
xcellList <- as.list(read.table("./data-raw/xCell_genes.txt",
                                sep = "\t",
                                header = T,
                                stringsAsFactors = F))

names(xcellList) <- sapply(strsplit(names(xcellList),"_"), function(x) x [[1]])

methodSignatures$xCell <- mergeDups(xcellList)

methodSignatures <- lapply(methodSignatures, function(geneSet) {
  lapply(geneSet, function(cellGenes) {
    cellGenes <- removeBlanks(cellGenes)
    cellGenes <- sort(cellGenes)
    return(cellGenes)
  })
})

# CIBERSORT LM22 Signature Matrix

methodSignatures$CIBERSORT <- read.table("./data-raw/CIBERSORT_LM22.txt",
                                 sep = "\t",
                                 header = T,
                                 row.names = "Gene.symbol")

# TIMER-Style Immune Filter, Gene Sets With Negative Purity Correlation  (All Cancers)

methodSignatures$ImmuneGenes <- read.table("./data-raw/Purity_immune_genes.txt",
                             sep = "\t",
                             header = T,
                             stringsAsFactors = F)

# ESTIMATE Gene Sets With Negative Purity Correlation (All Cancers)

methodSignatures$ESTIMATE <- read.table("./data-raw/ESTIMATE_stromal_genes.txt",
                                sep = "\t",
                                header = T,
                                stringsAsFactors = F)


save(methodSignatures, file = "./data/methodSignatures.rda")
