\name{CellCounts}
\alias{CellCounts}
\title{Compute Number of Cells Positive for Certain Cytokine Combinations}
\usage{
CellCounts(data, combinations)
}
\arguments{
  \item{data}{A list of matrices. Each matrix \code{i} is
  of dimension \code{N_i cells} by \code{K} common
  markers.}

  \item{combinations}{A list of 'combinations' == integer
  vectors, which is used to denote the subsets. See the
  examples.}
}
\description{
This function takes a list of cell expression matrices,
each matrix \code{i} of dimension \code{N_i cells} by
\code{K} common markers.
}
\examples{
K <- 6 ## number of markers
data <- replicate(10, simplify=FALSE, {
  m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
  m[m < 2500] <- 0
  colnames(m) <- c("IL2", "IL4", "IL6", "Mip1B", "IFNg", "TNFa")
  return(m)
})
names(data) <- sample(letters, 10)

combos <- list(1, 2, 3, 4, 5, 6) ## marginal cell counts
cc <- CellCounts(data, combos)
f <- function(data) {
  do.call(rbind, lapply(data, function(x) apply(x, 2, function(x) sum(x > 0))))
}
cc2 <- f(data)
stopifnot(identical( unname(cc), unname(cc2) ))

## We can also generate cell counts by expressing various combinations
## of markers (names) in the data.
CellCounts(data, "IL2|IL4") ## count cells expressing IL2 or IL4
CellCounts(data, "IL2|IL4|IL6") ## count cells expressing IL2, IL4 or IL6
CellCounts(data, c("IL2", "IL4", "IL6")) ## counts for each of IL2, IL4, IL6 (marginally)
}
\seealso{
\code{\link{discrete_combinations}}
}

