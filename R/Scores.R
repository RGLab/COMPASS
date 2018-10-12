##' Compute the Functionality Score for each subject fit in a COMPASS model
##'
##' Computes the functionality score for each observation from the gamma matrix
##' of a COMPASS model fit. The scores are normalized according to the total
##' number of possible subsets that could be observed (\code{2^M - 1}).
##'
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}. Alternatively, a matrix of functionality scores,
##'   used under the assumption that the 'null' category has been dropped.
##' @param n The number of markers included in an experiment. It is inferred
##'   from the data when \code{x} is a \code{COMPASSResult}.
##' @param markers The set of markers for which to compute a Functionality score. By default uses all markers. Must match names returned by \code{markers()}.
##' @return A numeric vector of functionality scores.
##' @export
##' @examples
##' FunctionalityScore(CR)
##' @note The null category is implicitly dropped when computing the functionality
##'   score for a \code{COMPASS} result; this is not true for the regular matrix
##'   method.
##' @import pdist
FunctionalityScore <- function(x, n,markers=NULL) {
  UseMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.COMPASSResult <- function(x, n,markers=NULL) {
  ## we drop the last column as it is the 'NULL' category
  n <- ncol(x$fit$categories) - 1
  y <- x$fit$mean_gamma[, -ncol(x$fit$mean_gamma),drop=FALSE]
  fs = apply(y, 1, function(row) {
    sum(row) / (2^n - 1)
  })
  #If markers was given, we compute a functionality score based on the subset of markers
  if(!is.null(markers)){
    if(!all(markers%in%markers(x))){
      stop("Invalid marker names")
    }
    message("Computing a Functionality Score based on ",paste(markers,collapse=", "))
    new_categories = unique(categories(x,FALSE)[,markers,drop=FALSE])
    all_categories=categories(x,FALSE)[,markers,drop=FALSE]
    dmat = as.matrix(pdist(new_categories,all_categories))
    cat_indices = apply(dmat,1,function(y)which(y==0))
    new_mean_gamma=apply(cat_indices,2,function(i)apply(Gamma(x)[,i,],1,mean))
    new_scores = rowSums(new_mean_gamma)
    new_scores=new_scores/(2^length(markers)-1)
    fs=new_scores
  }
  return(fs)
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.default <- function(x, n,markers=NULL) {
  apply(x, 1, function(row) {
    sum(row) / (2^n - 1)
  })
}

##' Compute the Polyfunctionality Score for each subject fit in a COMPASS model
##'
##' Computes the Polyfunctionality score for each observation from the
##' gamma matrix of a \code{COMPASS} model fit. The scores are normalized to
##' one.
##'
##' @param x An object of class \code{COMPASSResult}, as returned by
##'   \code{\link{COMPASS}}. Alternatively, a matrix of functionality scores.
##' @param degree A vector of weights. If missing, this is the 'degree of
##' functionality', that is, the number of markers expressed in a particular
##' category.
##' @param n The total number of markers. This is inferred when \code{x} is a
##'   \code{COMPASSResult}, and is unused in that case.
##' @param markers A \code{character} specifying the markers for which to compute the score. Must match names in \code{markers()}.
##' @return A numeric vector of polyfunctionality scores.
##' @export
##' @examples
##' PolyfunctionalityScore(CR)
PolyfunctionalityScore <- function(x, degree, n,markers=NULL) {
  UseMethod("PolyfunctionalityScore")
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.COMPASSResult <- function(x, degree, n,markers=NULL) {
  degree <- x$fit$categories[, "Counts"]
  n <- ncol(x$fit$categories) - 1
  y <- x$fit$mean_gamma
  pfs= apply(y, 1, function(row) {
    ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
    sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
  })
  if(!is.null(markers)){
    if(!all(markers%in%markers(x))){
      stop("Invalid marker names")
    }
    message("Computing a Polyfunctionality Score based on ",paste(markers,collapse=", "))
    new_categories = unique(categories(x,FALSE)[,markers,drop=FALSE])
    all_categories=categories(x,FALSE)[,markers,drop=FALSE]
    dmat = as.matrix(pdist(new_categories,all_categories))
    cat_indices = apply(dmat,1,function(y)which(y==0))
    new_mean_gamma=apply(cat_indices,2,function(i)apply(Gamma(x)[,i,],1,mean))
    degree <- rowSums(new_categories)
    n <- ncol(new_categories)
    y <- new_mean_gamma
    pfs= apply(y, 1, function(row) {
      ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
      sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
    })
  }
  return(pfs)
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.default <- function(x, degree, n,markers=NULL) {
  apply(x, 1, function(row) {
    ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
    sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
  })
}
