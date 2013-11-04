## debugging
if (FALSE) {
  library(flowWorkspace)
  gs <- load_gs("../shinyGate/gs/LoveLab-Gated")
}

##' Extract A flowWorkspace Node to a COMPASSContainer
##' 
##' This function can be used to take gated data from a \code{flowWorkspace}
##' object (e.g., a \code{GatingSet}), and retrieve it as a 
##' \code{COMPASSContainer}, suitable for ##' modelling with \code{COMPASS}.
##' 
##' @param gs A \code{GatingSet} object.
##' @param node The name of the node. We extract all children of that node,
##'   under the assumption that the markers of interest have been gated at
##'   the chosen node.
##' @param children A character vector, or list, giving the names of the children of
##'   \code{node} that we wish to extract. If unspecified, we take all
##'   children of that node. If the character vector passed is named, we map
##'   from value to name; e.g. if a node is named \code{154} but we prefer the
##'   name \code{CD154+}, we could pass \code{`CD154+`="154"}.
##' @param meta A \code{data.frame} of meta data. If it is missing, we attempt
##'   to pull the metadata from the \code{GatingSet}, and throw an error if it
##'   is unavailable.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn. 
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @export
##' @seealso \code{\link{COMPASSContainer}}
GatingSetToCOMPASS <- function(gs, node, children, 
  meta, individual_id, sample_id) {
  
  if (require(flowWorkspace)) {
    
    if (!inherits(gs, "GatingSet")) {
      stop("'gs' must be an object of class 'GatingSet'")
    }
    
    doWarn <- FALSE
    .children <- getChildren(gs[[1]], node, isPath=TRUE)
    for (i in 1:length(gs)) {
      if (!identical(.children, getChildren(gs[[i]], node, isPath=TRUE))) {
        doWarn <- TRUE
      }
    }
    
    if (doWarn) {
      warning("The children of node '", node, "' are not identical across ",
        "each of the samples encapsulated by your GatingSet; be sure that ",
        "your samples have been gated uniformly at each node.")
    }
    
    ## get the data at each node
    data <- vector("list", length(children))
    for (i in 1:length(gs)) {
      data[[i]] <- getData(gs, children[i])
    }
    data <- getData(gs, node)
    intensities <- Map(exprs, data)
    names(intensities) <- sampleNames(gs)
    counts <- integer( length(gs) )
    for (i in seq_along(counts)) {
      counts[[i]] <- getTotal( gs[[i]], node, flowJo=FALSE )
    }
    names(counts) <- sampleNames(gs)
    
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
    
    return(CC)
    
  } else {
    stop("This function requires 'flowWorkspace' to be installed; ",
      "try running 'biocLite(\"flowWorkspace\")' to install.")
  }
  
}
