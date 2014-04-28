## first, we generate the 'combinations' vector, which is in a specific
## format used by the C++ code 'CellCounts'; see 'CellCounts.cpp'
## for details
combos <- as.list( as.data.frame( apply(categories, 1, function(x) {
  tmp <- c( which(x == 1), -which(x == 0) )
  tmp <- tmp[ match(1:length(tmp), abs(tmp)) ]
  return(tmp)
})))

.counts <- function(y, combos, counts) {
  y <- lapply(y, function(yy) {
    apply(yy, 2, function(yyy) {
      yyy > 0
    })
  })
  m <- .Call("MIMOSA_CellCounts", y, combos, PACKAGE="MIMOSA")

  ## set the last column to be the 'null'
  m[, ncol(m)] <- counts[ names(y) ] - apply(m[,-ncol(m), drop=FALSE], 1, sum)
  return(m)
}

n_s <- .counts(y_s, combos, counts)
n_u <- .counts(y_u, combos, counts)
