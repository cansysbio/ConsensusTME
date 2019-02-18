#### Useful Functions ####

removeBlanks <- function(list.in) {
# Removes blanks from vector of character or a list containing gene sets
#  
# Args:
#   list.in: Either character vector or list containing gene sets
#
# Returns:
#   The vector or list with any blanks removed
  
  if (class(list.in) == 'list') {
    list.out <- lapply(list.in[list.in != ''], function(x) x[x != ''])
    return(list.out)
  } else if (class(list.in) == 'character') {
    list.out <- list.in[list.in != '']
    return(list.out)
  } else {
    sprintf('%s is not a class this function can use', class(list.in))
  }
}

mergeDups <- function(list.in, removenames = TRUE) {
# When supplied with a list with duplicate names, merges elements of list
# with an identical name.
# 
# Args:
#  list.in: List containing elements with the same name
#  removenames: if true removes the name of each string within list elements
# Returns:
#  A list with the elements with the same name combined

    if (class(list.in) == 'list'){
    lnames <- names(list.in)
    ulnames <- unique(lnames)
    list.grouped <- lapply(ulnames, function(x) list.in[lnames == x])
    list.merged <- lapply(list.grouped, unlist, recursive = F)
    if (removenames) {
      list.merged <- lapply(list.merged, unname)
    }
    names(list.merged) <- ulnames
    return (list.merged)
  } else {
    sprintf('%s is not a class this function can currently use',class(list.in))
  }
}

longToWide <- function(df.in, cellcol = NULL,  genecol = NULL) {
# From a dataframe with genes in one column and cell types in another create a
# dataframe with cell names as columns and genes underneath
#
# Args:
#  df.in: List containing two (or more columns), NB: Requires stringsAsFactors = T
#  cellcol: column name containing cell names (default 1st)
#  genecol: column name containing genes (default 2nd)
#
# Returns:
#  A dataframe containing the cell types as columns with genes beneath
  require(data.table)
  require(reshape2)
  
  if (ncol(df.in) > 2) {
    if (is.null(cellcol) | is.null(genecol)) {
      stop('More than two columns have been inputted without specifying
           which columns contain the genes and celltypes')
    } else if (is.numeric(c(cellcol, genecol))){
      cellcol = colnames(df.in)[cellcol]
      genecol = colnames(df.in)[genecol]
    }
    if (!(cellcol %in% colnames(df.in)) | !(genecol %in% colnames(df.in))){
      stop('Index or Name of Column Not in Inputted Dataframe')
    } 
    tmp.df <- df.in[ ,c(cellcol, genecol)]
  } else if (ncol(df.in) == 2) {
    tmp.df <- df.in
  } else {
    stop('Error with Input Dataframe')
  }
  
  colnames(tmp.df) <- c('cell','gene')
  tmp.df$ind <- with(tmp.df, ave(seq_along(cell), cell, FUN = seq_along))
  new.df <- dcast(tmp.df, ind~cell, value.var = 'gene', fill = '')[ ,-1]
  new.df[is.na(new.df)] = ''
  return(new.df)
}
