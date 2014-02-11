##' COMPASS (Combinatorial Polyfunctionality Analysis of Single-Cells)
##' 
##' This package implements a model for the analysis of polyfunctionality
##' in single-cell cytometry experiments. The model effectively
##' identifies combinations of markers that are differentially expressed
##' between samples of cells subjected to different stimulations.
##' 
##' @docType package
##' @name COMPASS-package
##' @useDynLib COMPASS, .registration=TRUE
##' @importFrom Rcpp evalCpp
##' @import data.table
##' @import grid
##' @alias COMPASS-package
##' @seealso
##'   \itemize{
##'   \item \code{\link{COMPASSContainer}}, for information on getting your
##'   cytometry data into a suitable format for use with \code{COMPASS},
##'   \item \code{\link{COMPASS}}, for the main model fitting routine.
##'   }
NULL
