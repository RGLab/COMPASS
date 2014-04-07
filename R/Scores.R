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
##' @return A numeric vector of functionality scores.
##' @export
##' @examples
##' FunctionalityScore(CR)
##' @note The null category is implicitly dropped when computing the functionality
##'   score for a \code{COMPASS} result; this is not true for the regular matrix
##'   method.
FunctionalityScore <- function(x, n) {
  UseMethod("FunctionalityScore")
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.COMPASSResult <- function(x, n) {
  ## we drop the last column as it is the 'NULL' category
  n <- ncol(x$fit$categories) - 1
  x <- x$fit$mean_gamma[, -ncol(x$fit$mean_gamma)]
  apply(x, 1, function(row) {
    sum(row) / (2^n - 1)
  })
}

##' @rdname FunctionalityScore
##' @export
FunctionalityScore.default <- function(x, n) {
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
##' @return A numeric vector of polyfunctionality scores.
##' @export
##' @examples
##' PolyfunctionalityScore(CR)
PolyfunctionalityScore <- function(x, degree, n) {
  UseMethod("PolyfunctionalityScore")
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.COMPASSResult <- function(x, degree, n) {
  degree <- x$fit$categories[, "Counts"]
  n <- ncol(x$fit$categories) - 1
  x <- x$fit$mean_gamma
  apply(x, 1, function(row) {
    ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
    sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
  })
}

##' @rdname PolyfunctionalityScore
##' @export
PolyfunctionalityScore.default <- function(x, degree, n) {
  apply(x, 1, function(row) {
    ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
    sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
  })
}
