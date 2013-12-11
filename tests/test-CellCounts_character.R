library(COMPASS)

n <- 1E2

data <- replicate(100, matrix( rnorm(n * 6, 2500, 500), nrow=n), simplify=FALSE)
data <- lapply(data, function(x) {
  colnames(x) <- LETTERS[1:6]
  x[ x < 1500 ] <- 0
  return (x)
})
combinations <- colnames(data[[1]])
CellCounts(data, combinations)
stopifnot( identical(
  unname(m1 <- CellCounts(data, combinations)),
  unname(m2 <- CellCounts(data, 1:6))
) )
