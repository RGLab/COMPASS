##' Extract A flowWorkspace Node to a COMPASSContainer
##' 
##' This function can be used to take gated data from a \code{flowWorkspace}
##' object (e.g., a \code{GatingSet}), and retrieve it as a 
##' \code{COMPASSContainer}, suitable for ##' modelling with \code{COMPASS}.
##' 
##' @param gs A \code{GatingSet} object.
##' @param node The name of the node.
##' @param meta A \code{data.frame} of meta data. If it is missing, we attempt
##'   to pull the metadata from the \code{GatingSet}, and throw an error if it
##'   is unavailable.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn. 
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @export
GatingSetToCOMPASS <- function(gs, node, meta, individual_id, sample_id) {
  
  if (require(flowWorkspace)) {
    
    if (!inherits(gs, "GatingSet")) {
      stop("'gs' must be an object of class 'GatingSet'")
    }
    
    data <- getData(gs, node)
    intensities <- Map(exprs, data)
    names(intensities) <- sampleNames(gs)
    counts <- unlist(lapply(gs, function(x) getTotal(x, "root", flowJo=FALSE)))
    
    ## try to get the metadata, if available
    if (missing(meta)) {
      tryCatch(meta <- pData(gs),
        error=function(e) {
          stop("No metadata available!")
        })
    }
    
    CC <- COMPASSContainer(
      data=intensities,
      counts=counts,
      meta=meta,
      individual_id=individual_id,
      sample_id=sample_id
    )
    
  } else {
    stop("This function requires 'flowWorkspace' to be installed; ",
      "try running 'biocLite(\"flowWorkspace\")' to install.")
  }
  
}
