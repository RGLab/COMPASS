library(testthat)
## functions common to one or more of 'preprocess.R', 'ui.R', 'server.R'

## get combinations up to an order n
## combine is a function that tells us how
## we combine a particular set of combos
## note that order of elements is 'preserved'
combinations <- function(x, combine=c, ...) {
  combine <- match.fun(combine)
  output <- vector("list", length(x))
  for( i in 1:length(x) ) {
    tmp <- combn(x, i, simplify=FALSE)
    output[[i]] <- lapply(tmp, combine, ...)
  }
  return( unlist(output, recursive=FALSE) )
}

stopifnot( identical( unlist( combinations(1:6)[1:6] ), 1:6 ) )

## for computing cytokine combinations
marginal <- function(dat, combos, ncol, colnames) {
  d <- vector("list", length(dat))
  for( i in 1:length(dat) ) {
    cat("Iteration", i, "of", length(dat), ".\n")
    x <- dat[[i]]
    output <- matrix(FALSE, nrow=nrow(x), ncol=ncol)
    for( j in 1:ncol(output) ) {
      output[, j] <- as.logical(apply( x[, combos[[j]], drop=FALSE], 1, prod))
    }
    colnames(output) <- colnames
    d[[i]] <- output
  }
  return(d)
}

joint <- function(dat, combos, ncol, colnames) {
  
  ## convert dat to binary
  dat[] <- lapply(dat, function(x) {
    colApply(x, as.logical, drop=FALSE)
  })
  
  d <- vector("list", length(dat))
  for( i in 1:length(dat) ) {
    cat("Iteration", i, "of", length(dat), ".\n")
    x <- dat[[i]]
    output <- matrix(FALSE, nrow=nrow(x), ncol=ncol)
    for( j in 1:ncol(output) ) {
      check <- rep(FALSE, ncol(x))
      check[ combos[[j]] ] <- TRUE
      output[, j] <- apply(x, 1, function(x) { all(x == check) })
    }
    colnames(output) <- colnames
    d[[i]] <- output
  }
  return(d)
}

## translates phenotype to label for different plots
phenoToLabel <- function(x) {
  return( switch(x,
    `MeanGamma`="Probability of Ag-Specificity",
    x
  ) )
}

## transform subset label from boolean-style to +/- style
transform_subset_label <- function(x) {
  x <- as.character(x)
  sapply(x, function(x) {
    splat <- unlist(strsplit(x, "&", fixed=TRUE))
    output <- lapply(splat, function(y) {
      if (substring(y, 1, 1) == "!") {
        return( paste0( substring(y, 2, nchar(y)), "-"))
      } else {
        return(paste0(y, "+"))
      }
    })
    do.call( function(...) paste(..., sep="", collapse=""), output )
  })
}
