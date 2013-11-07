##' Generate all Boolean Subsets
##' 
##' TODO
##' 
##' @param x A \code{COMPASSContainer}.
##' @param as.matrix Boolean; if \code{TRUE} we return results as a matrix;
##'   otherwise, we return the results as a list.
##' @export
BooleanSubsets <- function(x, as.matrix) {
  UseMethod("BooleanSubsets")
}

##' @rdname BooleanSubsets
##' @method BooleanSubsets COMPASSContainer
##' @S3method BooleanSubsets COMPASSContainer
BooleanSubsets.COMPASSContainer <- function(x, as.matrix=FALSE) {
  x <- x$data
  NextMethod("BooleanSubsets")
}

##' @rdname BooleanSubsets
##' @method BooleanSubsets default
##' @S3method BooleanSubsets default
BooleanSubsets.default <- function(x, as.matrix=FALSE) {
  
  combinations <- unique( as.data.table( 
    lapply( 
      lapply( 
        as.data.table( do.call( rbind, x ) ), 
        as.logical 
      ), 
      as.integer 
    ) 
  ) )
  
  combinations[, c("Counts") := apply(.SD, 1, sum)]
  setkeyv(combinations, c("Counts", rev(names(combinations))))
  combinations[, Counts := NULL]
  if (as.matrix) {
    combinations <- as.matrix(combinations)
  } else {
    combinations <- as.list( as.data.frame( apply(combinations, 1, function(x) {
      combinations <- c( which(x == 1), -which(x == 0) )
      combinations <- combinations[ match(1:length(combinations), abs(combinations)) ]
      return(combinations)
    })))
  }
  return(combinations)
}
