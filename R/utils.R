stopf <- function(...) stop(gettextf(...), call.=FALSE)
warnf <- function(...) warning(gettextf(...), call.=FALSE)

clamp <- function(x, from, to) {
  x[ x < from ] <- from
  x[ x > to ] <- to
  x
}

##' Markers
##'
##' Returns the markers associated with an experiment.
##'
##' @param object An \R object.
##' @export
markers <- function(object) {
  UseMethod("markers")
}

##' @export
markers.COMPASSContainer <- function(object) {
  colnames(object$data[[1]])
}

##' @export
markers.COMPASSResult <- function(object) {
  cats <- categories(object, FALSE)
  colnames(cats)
}

##' Categories
##'
##' Returns the categories matrix in a \code{COMPASSResult} object.
##'
##' @param x A \code{COMPASSResult} object.
##' @param counts Boolean; if \code{TRUE} we return the counts (degree of
##'   functionality) as well.
##' @export
categories <- function(x, counts) {
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' is not a 'COMPASSResult'", call.=FALSE)
  }
  cats <- x$data$categories
  if (!counts) {
    cats <- cats[, colnames(cats) != "Counts", drop=FALSE]
  }
  return(cats)
}

##' Metadata
##'
##' Functions for getting and setting the metadata associated with an object.
##'
##' @param x An \R object.
##' @param value An \R object appropriate for storing metadata in object \code{x};
##'   typically a \code{data.frame}.
##' @export
metadata <- function(x) UseMethod("metadata")

##' @export
metadata.COMPASSContainer <- function(x) {
  x$meta
}

##' @export
metadata.COMPASSResult <- function(x) {
  x$data$meta
}

##' @export
`metadata<-` <- function(x, value) UseMethod("metadata<-")

##' @export
`metadata<-.COMPASSContainer` <- function(x, value) {
  x$meta <- check_meta(x$data, x$counts, value, x$individual_id, x$sample_id)
  x
}

##' COMPASSResult Accessors
##'
##' These functions can be used for accessing data within a \code{COMPASSResult}.
##'
##' @param x A \code{COMPASSResult} object.
##' @name COMPASSResult-accessors
NULL

##' The gamma array associated with a \code{COMPASS} model fit.
##' @rdname COMPASSResult-accessors
Gamma <- function(x) {
  return(x$fit$gamma)
}

##' @rdname COMPASSResult-accessors
MeanGamma <- function(x) {
  return(x$fit$mean_gamma)
}

