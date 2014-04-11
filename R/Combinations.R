##' Generate Combinations
##' 
##' Given an intenger \code{n}, generate all binary combinations of
##' \code{n} elements, transformed to indices. This is primarily for use
##' with the \code{\link{CellCounts}} function, but may be useful for
##' users in some scenarios.
##' 
##' @param n An integer.
##' @export
##' @examples
##' Combinations(3)
Combinations <- function(n) {
  values <- do.call( function(...) {
    expand.grid(..., KEEP.OUT.ATTRS=FALSE)
  }, replicate(n, c(1, -1), simplify=FALSE))
  for (i in 1:ncol(values)) {
    values[, i] <- values[, i] * i
  }
  return( transpose_list(values) )
}
