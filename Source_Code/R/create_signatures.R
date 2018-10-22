createSignatures <- function(saveRdata = F, saveTxt = F, path.save.sigs = '.' ){
  # Creates a list containg cancer specific gene sets for each cell type, also
  # included in this list is an unspecific list which could be used for any analyses
  # on cancer types not included.
  #
  # Args:
  #  saveRdata: If True, the gene lists for all cancers will be saved as an RData file (in working directory by default)
  #  saveTxt: If True, the gene lists for all cancers will be saved as individual tab-seperated files
  #           within a directory named 'Individual-Cancer-Types', this will be created (in working directory be default)
  #  path.save.sigs: Allows user to specify a path for the gene signatures to be saved
  #
  # Returns:
  #   Returns a list containing cell-type specific gene lists for each cancer type along with a non-specific list of genes
  #   that has not been filtered using tumour purity

  source('~/Documents/ConsensusTME_Work/Consensus_Git/ConsensusTME/Source_Code/R/create_method_genes_object.R')

  #### Create Cancer Specific Consensus Gene List ####

  # Create List Of Genes For Each Method

  method.genes <- createMethodGenesList(saveRdata = F, path.save = NULL)

  # Create Consensus List Of Genes

  consensus.all <- vector('list', 18)
  cell.types <- unique(unlist(lapply(method.genes[!(names(method.genes) %in% c('TIMER','ESTIMATE'))], names)))
  names(consensus.all) <- cell.types

  for (cell in cell.types) {
    consensus.all[[cell]] <- unique(unlist(lapply(method.genes[!(names(method.genes) %in% c('TIMER','ESTIMATE'))]
                                                  ,function(x) x[names(x) == cell])))
  }

  consensus.all <- lapply(consensus.all, sort)

  # Filter Genes To Only Include Genes With Negative Purity Correlation For Each Cancer Type

  timer.all <- method.genes$TIMER
  estimate.all <- method.genes$ESTIMATE
  cancer.types <- unique(timer.all[ ,'Disease.Name'])
  cancer.types <- c(cancer.types, 'Unspecific')
  stromal.cells <- c('Endothelial','Fibroblasts')

  consensus.cancer <- setNames(vector('list',33), cancer.types)

  for (canc in cancer.types){
     immune.genes <- timer.all[timer.all$Disease.Name == canc, 1]
     stromal.genes <- estimate.all[estimate.all$Cancer == canc, 1]
     if (canc == 'Unspecific') {
       cancer.sig <- consensus.all
     } else {
       cancer.sig <- setNames(vector('list', 18), cell.types)
       for(cells in cell.types){
         if(cells %in% stromal.cells){
           cancer.sig[[cells]] <- intersect(consensus.all[[cells]], stromal.genes)
         } else {
           cancer.sig[[cells]] <- intersect(consensus.all[[cells]], immune.genes)
         }
       }
     }
     consensus.cancer[[canc]] <- cancer.sig
     if (saveTxt) {
       sig.df <- t(ldply(cancer.sig, rbind))
       sig.df[is.na(sig.df)] <- ''
       path.txt <- sprintf('%s/Individual_Cancer_Types', path.save.sigs)
       if (!dir.exists(path.txt)) {
         dir.create(path.txt)
       }
       write.table(sig.df, sprintf('%s/%s_Gene_Signatures.txt', path.txt, canc),
                   sep = '\t', row.names = F, col.names = F)
     }
  }

  # Save Gene Lists

  if (saveRdata) {
    save(consensus.cancer, file = sprintf('%s/ConsensusTME_GeneSets.RData', path.save.sigs))
  }
  return(consensus.cancer)
}
