##' Compute the Functionality Score for each subject fit in a COMPASS model
##' 
##' Computes the functionality score for each observation from the gamma matrix
##' of a COMPASS model fit.
##' 
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}. Alternatively, a matrix of functionality scores.
##' @return A numeric vector of functionality scores.
##' @export
FunctionalityScore <- function(x) {
  UseMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @method FunctionalityScore COMPASSResult
##' @S3method FunctionalityScore COMPASSResult
FunctionalityScore.COMPASSResult <- function(x) {
  x <- x$fit$mean_gamma
  NextMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @method FunctionalityScore default
##' @S3method FunctionalityScore default
FunctionalityScore.default <- function(x) {
  Fscore <- rowMeans(x[, -c(ncol(x))])
  return(Fscore)
}

##' Compute the Polyfunctionality Score for each subject fit in a COMPASS model
##' 
##' Computes the Polyfunctionality score for each observation from the 
##' gamma matrix of a COMPASS model fit.
##' 
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}.
##' @param normalization A \code{character} vector specifying how the score
##'  is to be normalized. Either using \code{"all"} possible categories, or 
##'  just the \code{"observed"} categories. Defaults to \code{"all"} categories.
##' @return a \code{vector} of polyfunctionality scores.
##' @export
PolyfunctionalityScore <- function(x, degree, normalization) {
  UseMethod("PolyfunctionalityScore")
}

##' @rdname PolyfunctionalityScore
##' @method PolyfunctionalityScore COMPASSResult
##' @S3method PolyfunctionalityScore COMPASSResult
PolyfunctionalityScore.COMPASSResult <- function(x, degree, normalization="all") {
  M <- x$fit$mean_gamma
  degree <- rev(rev(x$fit$categories[,ncol(x$fit$categories)])[-1L])
  switch(normalization,
    all = {
      norm <- choose(
        ncol(x$fit$categories)-1, 
        1:(ncol(x$fit$categories)-1)
      )
    },
    observed = {
      norm <- table(degree)[degree]
    }
  )
  degree <- degree / norm[degree]
  PFscore <- M[,-ncol(M)] %*% degree
  PFscore <- setNames( as.numeric(PFscore), rownames(PFscore) )
  return(PFscore)
}

##' @rdname PolyfunctionalityScore
##' @method PolyfunctionalityScore COMPASSResult
##' @S3method PolyfunctionalityScore COMPASSResult
PolyfunctionalityScore.default <- function(x, degree, normalization) {
  return( rowMeans(as.matrix(x) * matrix( rep(degree, each=nrow(x)), nrow=nrow(x))) )
}
