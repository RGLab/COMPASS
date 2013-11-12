##' Subset a COMPASSContainer
##' 
##' Use this function to subset a \code{COMPASSContainer}.
##' 
##' @param x A \code{COMPASSContainer}.
##' @param subset A logical expression, evaluated within the metadata,
##'   indicating which entries to keep.
##' @param ... Optional arguents; currently unused.
##' @method subset COMPASSContainer
##' @S3method subset COMPASSContainer
subset.COMPASSContainer <- function(x, subset, ...) {
  
  call <- match.call()$subset
  
  keep <- eval(call, envir=x$meta)
  if (!is.logical(keep)) {
    stop("'subset' should evaluate to a logical result")
  }
  
  good_samples <- unique(x$meta[[ x$sample_id ]][keep])
  
  ## subset everything
  output <- COMPASSContainer(
    data=x$data[good_samples],
    counts=x$counts[good_samples],
    meta=x$meta[keep, ],
    individual_id=x$individual_id,
    sample_id=x$sample_id,
    stimulation_id=x$stimulation_id
  )
  setattr(output, "class", "COMPASSContainer")
  return(output)
  
}

