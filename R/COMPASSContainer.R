##' Generate the Data Object used by COMPASS
##'
##' This function generates the data container suitable for use with
##' \code{COMPASS}.
##'
##' The \code{names} attributes for the \code{data} and \code{counts}
##' objects passed should match.
##'
##' @param data A list of matrices. Each matrix \code{M_i} is made up of
##'   \code{N_i} cells by \code{K} markers; for example, it could be the
##'   intensity information from an intracellular cytokine experiment.
##'   Each element of the list should be named; this name denotes which
##'   sample the cell intensities were measured from.
##' @param counts A named integer vector of the cell counts(of the parent population) for each
##'   sample in \code{data}.
##' @param meta A \code{data.frame} of metadata, describing the individuals
##'   in the experiment. Each row in \code{meta} should correspond to a row
##'   in \code{data}. There should be one row for each sample;
##'   i.e., one row for each element of \code{data}.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn.
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @param countFilterThreshold Numeric; if the number of cells expressing at
##'   least one marker of interest is less than this threshold, we remove that
##'   file. Default is 0, which means filter is disabled.
##' @return A \code{COMPASSContainer} returns a list made up of the same
##' components as input the model, but checks and sanitizes the supplied data
##' to ensure that it conforms to the expectations outlied above.
##'
##' @export
##' @example examples/GenerateSampleCOMPASSContainer.R
COMPASSContainer <- function(data, counts, meta,
                             individual_id, sample_id, countFilterThreshold = 0) {

  ## check names
  if (is.null(names(data)))
    stop( "'", deparse(substitute(data)), "' must have its 'names' ",
          "attribute set", call.=FALSE)

  if (is.null(names(counts)))
    stop("'", deparse(substitute(counts)), "' must have its 'names' ",
         "attribute set", call.=FALSE)

  if (is.list(counts)) {
    counts <- unlist(counts)
  }
  storage.mode(counts) <- "integer"

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

  ## guess the number of markers
  ## we have to do some ugly stuff here since there's no guarantee that
  ## the user has actually passed in a list of matrices (could be
  ## NULL, NA, a named vector...) but we want to be helpful
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
  if (any(null_data)) {
    for (i in which(null_data)) {
      data[[i]] <- empty_matrix
    }
  }

  ## ensure that the number of markers, names are identical
  n_markers_all <- sapply(data, ncol)
  if (length(unique(n_markers_all)) != 1) {
    stop("The number of markers is not identical across samples.", call.=FALSE)
  }

  names_all <- lapply(data, colnames)
  if (!identical( Reduce(union, names_all), Reduce(intersect, names_all))) {
    stop("The marker names are not identical across samples.", call.=FALSE)
  }

  ## type checking
  .type_check <- function(y) {
    yy <- deparse(substitute(y))
    if (is.null(names(y)))
      stop("'", yy, "' must have the 'names' attribute set", call.=FALSE)
    if (!is.list(y))
      stop("'", yy, "' must be a named list", call.=FALSE)
    if (!all(sapply(y, is.matrix)))
      stop("All elements in '", yy, "' must be matrices", call.=FALSE)
    if (any(sapply(y, function(x) is.null(colnames(x)))))
      stop("All matrices in '", yy, "' must have column names", call.=FALSE)
    nm <- colnames(y[[1]])
    for (i in seq_along(y)) {
      if (!all(nm %in% colnames(y[[i]])))
        stop("All matrices in '", yy,
             "' must share the same column (marker) names", call.=FALSE)
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
      stop("Not all names in 'counts' are in the names of '", yy, "'.",
           call.=FALSE)
    }
    return( invisible(NULL) )
  }

  .check_has_names(data, counts)
  
  if(countFilterThreshold > 0){
    message("Filtering low counts")
    filter <- counts > countFilterThreshold
    keep.names <- names(counts)[filter]
    data <- data[keep.names]
    counts <- counts[keep.names]
    meta <- subset(meta, eval(as.name(sample_id)) %in% keep.names)
    message(gettextf("Filtering %s samples due to low counts", length(filter) -
                length(keep.names)))  
  }
  
  
  ## ensure that the counts are >= the number of rows in the data
  if (any(sapply(data, nrow) > counts)) {
    stop("There are entries in 'counts' that are greater than the ",
         "number of rows included in the 'data' matrices.", call.=FALSE)
  }

  ## check for negative values in the data
  if (any(sapply(data, function(x) any(x < 0)))) {
    warning("There appear to be negative intensities in the 'data' supplied.")
  }

  ## check the metadata
  meta <- check_meta(data, counts, meta, individual_id, sample_id)

  output <- list(
    data=data,
    counts=counts,
    meta=meta,
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
##' @export
##' @examples
##' print(CC)
print.COMPASSContainer <- function(x, ...) {
  cat("A COMPASSContainer with ",
      length(x$data), " samples from ", length(unique(x$meta[[ x$individual_id ]])),
      " individuals, containing data across ",
      ncol(x$data[[1]]), " markers.\n", sep='', ...)
}

##' Summarize a COMPASSContainer Object
##'
##' This function prints summary information about a \code{COMPASSContainer}
##' object -- the number of samples, basic information about the metadata,
##' and so on.
##'
##' @param object An object of class \code{COMPASSContainer}.
##' @param ... Optional arguments; currently ignored.
##' @export
##' @examples
##' summary(CC)
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
