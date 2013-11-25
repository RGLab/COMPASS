##' Generate the Data Object used by COMPASS
##' 
##' This function generates the data container suitable for use with
##' \code{COMPASS}.
##' 
##' @param data A list of matrices. Each matrix \code{M_i} is made up of 
##'   \code{N_i} cells by \code{K} markers; for example, it could be the 
##'   intensity information from an intracellular cytokine experiment. 
##'   Each element of the list should be named; this name denotes which 
##'   sample the cell intensities were measured from.
##' @param counts A named integer vector of the cell counts for each
##'   sample in \code{data}. This element is required if 'null'
##'   cells; i.e., cells that expressed no markers, are not included in
##'   \code{data}.
##' @param meta A \code{data.frame} of metadata, describing the individuals
##'   in the experiment. Each row in \code{meta} should correspond to a row
##'   in \code{data}. There should be one row for each sample;
##'   i.e., one row for each element of \code{data}.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn. 
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @param stimulation_id The name of the vector in \code{meta} that denotes
##'   the type of stimulation each sample received.
##' @export
COMPASSContainer <- function(data, counts, meta, 
  individual_id, sample_id, stimulation_id) {
  
  ## check names
  if (is.null(names(data)))
    stop( "'", deparse(substitute(data)), "' must have its 'names' ",
      "attribute set")
  
  if (is.null(names(counts)))
    stop("'", deparse(substitute(counts)), "' must have its 'names' ",
      "attribute set")
  
  ## remove NULLs
  null_data <- names(data)[sapply(data, is.null)]
  if (length(null_data)) {
    warning("The following samples had no cytometry information available:\n\t",
      paste(null_data, collapse=", "))
  }
  
  ## convert named vectors to matrices
  for (i in seq_along(data)) {
    if (!is.null(names(data[[i]]))) {
      nm <- names( data[[i]] )
      data[[i]] <- matrix( data[[i]], nrow=1 )
      colnames( data[[i]] ) <- nm
    }
  }
  
  ## convert NULLs to 0-row matrices
  n_markers <- NULL
  marker_names <- NULL
  for (i in 1:length(data)) {
    if (is.matrix(data[[i]])) {
      n_markers <- ncol(data[[i]])
      marker_names <- colnames(data[[i]])
      break
    }
  }
  
  ## replace null with 0-row matrices
  empty_matrix <- matrix(nrow=0, ncol=n_markers)
  colnames(empty_matrix) <- marker_names
  mode(empty_matrix) <- "numeric"
  null_data <- sapply(data, is.null)
  data[ sapply(data, is.null) ] <- empty_matrix
  
  ## type checking
  .type_check <- function(y) {
    yy <- deparse(substitute(y))
    if (is.null(names(y))) stop("'", yy, "' must have the 'names' attribute set")
    if (!is.list(y)) stop("'", yy, "' must be a named list")
    if (!all(sapply(y, is.matrix)))
      stop("All elements in '", yy, "' must be matrices")
    if (any(sapply(y, function(x) is.null(colnames(x)))))
      stop("All matrices in '", yy, "' must have column names")
    nm <- colnames(y[[1]])
    for (i in seq_along(y)) {
      if (!all(nm %in% colnames(y[[i]])))
        stop("All matrices in '", yy, "' must share the same column (marker) names")
    }
  }
  
  .type_check(data)
  
  ## reorder y_s, y_u, so that the column names are in the same order
  .reorder_columns <- function(y) {
    nm <- colnames(y[[1]])
    for (i in seq_along(y)[-1]) {
      y[[i]] <- y[[i]][, match(colnames(y[[i]]), nm), drop=FALSE]
    }
    return(y)
  }
  
  data <- .reorder_columns(data)
  
  ## TODO: allow counts to be empty / NULL?
  
  ## ensure that all the names in counts are also in y_s, y_u
  .check_has_names <- function(y, counts) {
    yy <- deparse(substitute(y))
    if (!all(names(y) %in% names(counts))) {
      stop("Not all names in 'counts' are in the names of '", yy, "'.")
    }
    return( invisible(NULL) )
  }
  
  .check_has_names(data, counts)
  
  ## check the metadata
    
  if (!is.data.frame(meta))
    stop("'meta' must be a 'data.frame'")
  
  ## ensure that the names in y_s, y_u are also in the names of meta
  all_names <- names(data)
  missing_names <- all_names[ !(all_names %in% meta[[sample_id]]) ]
  
  if (length(missing_names)) {
    warning("There are ids in 'data' that are not available in ",
      "the metadata passed. These ids are: ", paste(missing_names, collapse=", "))
  }
  
  ## convert the sample id, individual id to character if necessary
  meta[[individual_id]] <- as.character( meta[[individual_id]] )
  meta[[sample_id]] <- as.character( meta[[sample_id]] )
  
  output <- list(
    data=data,
    counts=counts,
    meta=meta,
    stimulation_id=stimulation_id,
    individual_id=individual_id,
    sample_id=sample_id
  )
  
  class(output) <- "COMPASSContainer"
  return(output)
  
}

##' Print a COMPASSContainer Object
##' 
##' This function prints a \code{COMPASSContainer} object, giving basic
##' information about the object and the data it encapsulates.
##' 
##' @param x An object of class \code{COMPASSContainer}.
##' @param ... Optional arguments passed to \code{cat}.
##' @method print COMPASSContainer
##' @S3method print COMPASSContainer
##' @export
print.COMPASSContainer <- function(x, ...) {
  cat("An object of class COMPASSContainer with ", 
    length(x$data), " samples and ", 
    ncol(x$data[[1]]), " markers.\n", sep='', ...)
}

##' Summarize a COMPASSContainer Object
##' 
##' This function prints summary information about a \code{COMPASSContainer}
##' object -- the number of samples, basic information about the metadata,
##' and so on.
##' 
##' @method summary COMPASSContainer
##' @S3method summary COMPASSContainer
##' @param object An object of class \code{COMPASSContainer}.
##' @param ... Optional arguments; currently ignored.
##' @export
summary.COMPASSContainer <- function(object, ...) {
  print(object)
  cat( sep='', "The metadata describes ", ncol(object$meta), " variables on ",
    nrow(object$meta), " samples. The columns are:\n\n\t", 
    paste(names(object$meta), collapse=" "), "\n")
  hist( sapply(object$data, length),
    xlab="Number of Cells",
    ylab="Count",
    main="Number of Cells by Sample"
  )
}
