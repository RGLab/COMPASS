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
##' @export
##' @import plyr
##' @import data.table
##' @examples
##' scores(CR)
scores = function(x){
  if(class(x)!="COMPASSResult"){
    stop("x must be of class COMPASSResult")
  }
	FS = ldply(FunctionalityScore(x))
	PFS = ldply(PolyfunctionalityScore(x))
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
