\name{discrete_combinations}
\alias{discrete_combinations}
\title{Generate all possible combinations of integers from -n:n}
\usage{
discrete_combinations(n)
}
\arguments{
  \item{n}{An integer.}
}
\description{
This is primarily for use with the \code{\link{CellCounts}}
function.
}
\examples{
## generate some example data
K <- 6 ## number of markers
data <- replicate(10, simplify=FALSE, {
  m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
  m[m < 2500] <- 0
  colnames(m) <- c("IL2", "IL4", "IL6", "Mip1B", "IFNg", "TNFa")
  return(m)
})
names(data) <- sample(letters, 10)
str(data)
CellCounts(data, discrete_combinations(K)) ## all possible combos
}

