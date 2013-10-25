##' Compute the Functionality Score for each subject fit in a COMPASS model
##' 
##' Computes the functionality score for each observation from the gamma matrix
##' of a COMPASS model fit.
##' 
##' @param CR An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}.
##' @return A numeric vector of functionality scores.
##' @export
FunctionalityScore <- function(CR) {
  M <- CR$fit$mean_gamma
  Fscore <- rowMeans(M[, -c(ncol(M))])
  return(Fscore)
}

##' Compute the Polyfunctionality Score for each subject fit in a COMPASS model
##' 
##' Computes the Polyfunctionality score for each observation from the 
##' gamma matrix of a COMPASS model fit.
##' 
##' @param CR An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}.
##' @param normalization A \code{character} vector specifying how the score
##'  is to be normalized. Either using \code{"all"} possible categories, or 
##'  just the \code{"observed"} categories. Defaults to \code{"all"} categories.
##' @return a \code{vector} of polyfunctionality scores.
##' @export
PolyfunctionalityScore <- function(CR, normalization=c("all","observed")) {
  M <- CR$fit$mean_gamma
  degree <- rev(rev(CR$fit$categories[,ncol(CR$fit$categories)])[-1L])
  normalization <- match.arg(arg=normalization,choices=c("all","observed"))
  switch(normalization,
    all = {
      norm <- choose(
        ncol(CR$fit$categories)-1, 
        1:(ncol(CR$fit$categories)-1)
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
