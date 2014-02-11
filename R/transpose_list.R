##' Transpose a List
##' 
##' Transpose a matrix-like list.
##' @param x An \R list.
##' @export
##' @examples
##' l <- list( 1:3, 4:6, 7:9 )
##' stopifnot( identical(
##'   transpose_list( transpose_list(l) ), l
##' ) )
transpose_list <- function(x) {
  return( .Call(C_transpose_list, as.list(x)) )
}
