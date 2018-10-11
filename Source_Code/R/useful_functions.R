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

thisFile <- function() {
# Function from SO that allows a path for the script to reference itself
#
# Args: N/A
#
# Returns:
#  The path that the function script is in allowing reference to other datasets 

  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}