##' Generate Boolean Subsets
##' 
##' See the S3 dispatch methods for more information.
##' 
##' \itemize{
##'   \item{\code{\link{BooleanSubsets.COMPASSContainer}}}{ for a \code{COMPASSContainer}},
##'   \item{\code{\link{BooleanSubsets.default}}}{ for the default method}
##' }
##' 
##' @param x A \code{COMPASSContainer}.
##' @param as.matrix Boolean; if \code{TRUE} we return results as a matrix;
##'   otherwise, we return the results as a list.
##' @export
BooleanSubsets <- function(x, as.matrix) {
  UseMethod("BooleanSubsets")
}

##' Generate all available Boolean Subsets for a COMPASSContainer
##' 
##' Given a \code{COMPASSContainer}, we can generate all boolean
##' combinations of markers available in the data.
##' 
##' @rdname BooleanSubsets-COMPASSContainer
##' @method BooleanSubsets COMPASSContainer
##' @S3method BooleanSubsets COMPASSContainer
BooleanSubsets.COMPASSContainer <- function(x, as.matrix=FALSE) {
  x <- x$data
  NextMethod("BooleanSubsets")
}

##' Generate all available Boolean Subsets
##' 
##' Given a list of matrices, each with the same number of columns,
##' we generate all possible boolean subsets present in the data.
##' 
##' @rdname BooleanSubsets-default
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
