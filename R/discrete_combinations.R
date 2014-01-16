##' Generate all possible combinations of integers from -n:n
##' 
##' This is primarily for use with the \code{\link{CellCounts}} function.
##' 
##' @param n An integer.
##' @export
##' @examples
##' ## generate some example data
##' K <- 6 ## number of markers
##' data <- replicate(10, simplify=FALSE, {
##'   m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
##'   m[m < 2500] <- 0
##'   colnames(m) <- c("IL2", "IL4", "IL6", "Mip1B", "IFNg", "TNFa")
##'   return(m)
##' })
##' names(data) <- sample(letters, 10)
##' str(data)
##' CellCounts(data, discrete_combinations(K)) ## all possible combos
discrete_combinations <- function(n) {
  values <- do.call(expand.grid, replicate(n, c(1, -1), simplify=FALSE))
  for (i in 1:ncol(values)) {
    values[, i] <- values[, i] * i
  }
  return( unname(as.list( as.data.frame(t(values)) )) )
}
