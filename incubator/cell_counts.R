.check_cytokine_list <- function(data) {
  ncol <- sapply(data, ncol)
  stopifnot( length(unique(ncol)) == 1 )
  
  nm <- colnames(data[[1]])
  for (i in seq_along(data))
    if (!all(nm %in% colnames(data[[i]])))
      stop("All matrices must share the same column (marker) names")
  
  data <- lapply(data, function(x) {
    x[, match(colnames(x), nm)]
  })
  
  return(data)
  
}

setGeneric("cytokine_counts", function(data, combinations, ...) {
  standardGeneric("cytokine_counts")
})



setMethod("cytokine_counts", c("list", "character"), function(data, combinations, ...) {
  
  data <- .check_cytokine_list(data)
  c_names <- unlist( strsplit(combinations, "[-\\+]") )
  
})

##' Compute Number of Cells Positive for Certain Cytokine Combinations
##' 
##' This function takes a list of cell expression matrices, each matrix \code{i} 
##' of dimension \code{N_i cells} by \code{K} common markers.
##' 
##' @param data A list of matrices. Each matrix \code{i} is
##'   of dimension \code{N_i cells} by \code{K} common markers.
##' @param combinations A vector of