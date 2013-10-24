#'Compute the Functionality Score for each subject fit in a COMPASS model
#'
#'Computes the Functionality score for each observation from the gamma matrix of a COMPASS model fit.
#'
#'@param fit is a \code{COMPASSFit} object from either the discrete or continuous COMPASS model.
#'@return a \code{vector} of functionality scores.
#'@export
FunctionalityScore <- function(fit){
  M <- fit$mean_gamma
  Fscore <- rowMeans(M[,-c(ncol(M))])
  return(Fscore)
}

#'Compute the Polyfunctionality Score for each subject fit in a COMPASS model
#'
#'Computes the Polyfunctionality score for each observation from the gamma matrix of a COMPASS model fit.
#'
#'@param fit is a \code{COMPASSFit} object from either the discrete or continuous COMPASS model.
#'@param normalization a \code{character} vector specifying how the score is to be normalized. Either using \code{"all"} possible categories, or just the \code{"observed"} categories.
#'Defaults to \code{"all"} categories.
#'@return a \code{vector} of polyfunctionality scores.
#'@export
PolyfunctionalityScore <- function(fit,normalization=c("all","observed")){
  M <- fit$mean_gamma
  degree <- rev(rev(fit$categories[,ncol(fit$categories)])[-1L])
  normalization<-match.arg(arg=normalization,choices=c("all","observed"))
  switch(normalization,
  all = {norm <- choose(ncol(fit$categories)-1,1:(ncol(fit$categories)-1))
         },
  observed = { norm <- table(degree)[degree]
          })
  degree <- degree/norm[degree]
  PFscore <- M[,-ncol(M)]%*%degree
  return(PFscore)
}