## generate discrete combinations, as used for the CellCounts
## function
discrete_combinations <- function(n) {
  values <- do.call(expand.grid, replicate(n, c(1, -1), simplify=FALSE))
  for (i in 1:ncol(values)) {
    values[, i] <- values[, i] * i
  }
  return( unname(as.list( as.data.frame(t(values)) )) )
}
