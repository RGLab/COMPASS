##' Merge Two COMPASSContainers
##' 
##' This function merges two \code{COMPASSContainers}.
##' 
##' @param x A \code{COMPASSContainer}.
##' @param y A \code{COMPASSContainer}.
##' @param ... Optional arguments, currently ignored.
##' @export
##' @importFrom plyr rbind.fill
##' @examples
##' ## Chop the example COMPASSContainer into two, then merge it back together
##' CC1 <- subset(CC, trt == "Control")
##' CC2 <- subset(CC, trt == "Treatment")
##' merged <- merge(CC1, CC2)
##' res <- identical(CC, merge(CC1, CC2)) ## should return TRUE in this case
##' stopifnot( isTRUE(res) )
merge.COMPASSContainer <- function(x, y, ...) {
  
  nx <- names(x$data)
  ny <- names(y$data)
  
  if (length(ix <- intersect(nx, ny))) {
    stop("The sample names between the two COMPASSContainers are not unique!")
  }
  
  mx <- colnames(x$data[[1]])
  my <- colnames(y$data[[1]])
  
  ict <- intersect(mx, my)
  unn <- union(mx, my)
  
  if (any(!(unn %in% ict))) {
    warning("There are markers not in common between the two ",
      "COMPASSContainers and they will be dropped: ",
      paste(unn[ !(unn %in% ict) ], collapse=", ")
    )
  }
  
  ## Ensure the markers are in the same order
  x$data <- lapply(x$data, function(.) {
    return( .[, ict, drop=FALSE] )
  })
  
  y$data <- lapply(y$data, function(.) {
    return( .[, ict, drop=FALSE] )
  })
  
  data <- c(x$data, y$data)
  
  ## drop rows that are all 0s, if necessary
  data <- lapply(data, function(x) {
    x[ rowSums(x) != 0, , drop=FALSE ]
  })
  
  ## Merge the metadata
  meta <- rbind.fill(x$meta, y$meta)
  
  ## Merge the counts
  counts <- c(x$counts, y$counts)
  
  ## Check the metadata information
  stopifnot(x$individual_id == y$individual_id)
  stopifnot(x$sample_id == y$sample_id)
  
  ## Make a new COMPASSContainer
  return( COMPASSContainer(
    data=data,
    counts=counts,
    meta=meta,
    individual_id=x$individual_id,
    sample_id=x$sample_id
  ) )
  
}
