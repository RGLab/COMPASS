##' Print a COMPASSResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{COMPASS}} call.
##'
##' @param x An object of class \code{COMPASSResult}.
##' @param ... Optional arguments; currently unused.
##' @export
##' @examples
##' print(CR)
print.COMPASSResult <- function(x, ...) {
  n <- nrow(x$data$n_s)
  cat("A COMPASS model fit on", n, "paired samples.\n")
}

##' Fetch the table of scores and metadata from a COMPASSResult Object
##'
##' This function extracts the functionality and polyfunctionality scores from
##' a COMPASS result merged with the sample metadata table, accounting for any
##' dropped samples due to filtering.
##' @param x A \code{COMPASSResult} object.
##' @param markers A \code{character} vector of markers for which to compute the scores. Defaults to all markers. Must match the names returned by \code{markers()}.
##' @export
##' @import data.table
##' @examples
##' scores(CR)
scores = function(x,markers=NULL){
  if(class(x)!="COMPASSResult"){
    stop("x must be of class COMPASSResult")
  }
	FS = ldply(FunctionalityScore(x,markers=markers))
	PFS = ldply(PolyfunctionalityScore(x,markers=markers))
	m = x$data$meta
	setDT(m)
	setDT(FS)
	setDT(PFS)
	setnames(FS,c(".id","V1"),c(x$data$individual_id,"FS"))
	setnames(PFS,c(".id","V1"),c(x$data$individual_id,"PFS"))
	merge(merge(FS,PFS,by=x$data$individual_id),m,by=x$data$individual_id)
}

##' Summarize a COMPASSResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{COMPASS}} call.
##'
##' @param object An object of class \code{COMPASSResult}.
##' @param ... Optional arguments; currently unused.
##' @export
##' @examples
##' print(CR)
summary.COMPASSResult <- function(object, ...) {
  n <- nrow(object$data$n_s)
  cat("A COMPASSResult model fit on", n, "paired samples.\n")
}

#' Get a data.table of counts of polyfunctional subsets
#'
#'@param object An object of class \code{COMPASSResult}
#'@export
#'@import data.table methods
#'@importFrom grDevices bmp colorRamp dev.off hsv jpeg pdf png rainbow rgb rgb2hsv tiff
#'@importFrom graphics hist strwidth
#'@importFrom methods is
#'@importFrom stats as.dist cor dist hclust kmeans median na.omit runif sd setNames
#'@importFrom utils combn head installed.packages sessionInfo
#'@examples
#'getCounts(CR)
getCounts <- function(object){
  if(class(object)!="COMPASSResult"){
    stop("object must be of class COMPASSResult")
  }
  s = reshape2::melt(object$data$n_s)
  u = reshape2::melt(object$data$n_u)
  colnames(s) = c("ID","subset","stim")
  colnames(u) = c("ID","subset","unstim")
  counts = merge(s,u,by=c("ID","subset"))
  counts = data.table(counts)
  return(counts)
}
