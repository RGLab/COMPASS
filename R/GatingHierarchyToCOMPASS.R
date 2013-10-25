##' Extract A flowWorkspace Node to a COMPASSContainer
##' 
##' This function can be used to take gated data from a \code{flowWorkspace}
##' object (e.g., a \code{GatingHierarchy}), and retrieve it as a 
##' \code{COMPASSContainer}, suitable for ##' modelling with \code{COMPASS}.
##' 
##' @param gh A \code{GatingHierarchy} object.
##' @param node The name of the node.
##' @export
GatingHierarchyToCOMPASS <- function(gh, node) {
  if (!inherits(gh, "GatingHierarchy")) {
    stop("'gh' must be an object of class 'GatingHierarchy'")
  }
  
  data <- getData(gh, node)
  counts <- lapply(gh, function(x) getPopStats(x)[[3]])
  rownames(counts) <- rownames( getPopStats(gh[[1]]) )
  
}