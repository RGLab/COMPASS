##' Compute the Functionality Score for each subject fit in a COMPASS model
##' 
##' Computes the functionality score for each observation from the gamma matrix
##' of a COMPASS model fit.
##' 
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}. Alternatively, a matrix of functionality scores,
##'   used under the assumption that the 'null' category has been dropped.
##' @return A numeric vector of functionality scores.
##' @export
##' @examples
##' FunctionalityScore(CR)
FunctionalityScore <- function(x) {
  UseMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.COMPASSResult <- function(x) {
  ## we drop the last column as it is the 'NULL' category
  x <- x$fit$mean_gamma[, -ncol(x$fit$mean_gamma)]
  NextMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.default <- function(x) {
  return( rowMeans(x) )
}

##' Compute the Polyfunctionality Score for each subject fit in a COMPASS model
##' 
##' Computes the Polyfunctionality score for each observation from the 
##' gamma matrix of a \code{COMPASS} model fit.
##' 
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}. Alternatively, a matrix of functionality scores.
##' @param degree A vector of weights. If missing, this is the 'degree of
##' functionality', that is, the number of markers expressed in a particular
##' category.
##' @param normalization A \code{character} vector specifying how the score
##'  is to be normalized. Either using \code{"all"} possible categories, or 
##'  just the \code{"observed"} categories. Defaults to \code{"all"} categories.
##' @return A numeric vector of polyfunctionality scores.
##' @export
##' @examples
##' PolyfunctionalityScore(CR)
PolyfunctionalityScore <- function(x, degree, normalization) {
  UseMethod("PolyfunctionalityScore")
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.COMPASSResult <- function(x, degree, normalization="all") {
  M <- x$fit$mean_gamma
  if (missing(degree)) {
    degree <- x$fit$categories[,ncol(x$fit$categories)]
    degree <- degree[ -length(degree) ] ## drop the NULL category
  }
  switch(normalization,
    all = {
      norm <- choose(
        ncol(x$fit$categories)-1, 
        1:(ncol(x$fit$categories)-1)
      )
    },
    observed = {
      norm <- unname(c(table(degree)[degree]))
    }
  )
  degree <- degree / norm[degree]
  PFscore <- M[, -ncol(M)] %*% degree
  PFscore <- setNames( as.numeric(PFscore), rownames(PFscore) )
  
  ## divide by number of markers
  PFscore <- PFscore / (ncol(x$data$categories)-1)
  return(PFscore)
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.default <- function(x, degree, normalization) {
  return( rowMeans(as.matrix(x) * matrix( rep(degree, each=nrow(x)), nrow=nrow(x))) )
}
