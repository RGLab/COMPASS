##' Compute Total Cell Counts
##' 
##' This function is used to compute total cell counts, per individual,
##' from a \code{COMPASSContainer}.
##' 
##' @param data A \code{COMPASSContainer}.
##' @param subset An expression, evaluated within the metadata, defining
##'   the subset of \code{data} over which the counts are computed. If left
##'   unspecified, the counts are computed over all samples.
##' @param aggregate Boolean; if \code{TRUE} we sum over the individual,
##'   to get total counts across samples for each individual.
##' @export
##' @examples
##' TotalCellCounts(CC, trt == "Treatment")
##' TotalCellCounts(CC, trt == "Control")
##' TotalCellCounts(CC)
TotalCellCounts <- function(data, subset, aggregate=TRUE) {
  
  ## R CMD check silencers
  Counts <- PTID <- NULL
  
  if (!inherits(data, "COMPASSContainer"))
    stop("'data' must be an object of class 'COMPASSContainer'")
  
  subset_call <- match.call()$subset
  
  if (is.null(subset_call)) { ## implies subset is missing
    keep <- rep(TRUE, nrow(data$meta))
  } else {
    keep <- eval(subset_call, envir=data$meta)
  }
  
  samples_keep <- unique(data$meta[[ data$sample_id ]][ keep ])
  dt <- data.table(
    Samples=samples_keep
  )
  
  ## Merge in the counts
  counts_dt <- data.table(
    Samples=names(data$counts),
    Counts=unname(data$counts)
  )
  
  dt <- merge(dt, counts_dt, all.x=TRUE, by="Samples")
  
  if (aggregate) {
    
    ## Merge in the PTID information
    ptid_dt <- data.table(
      Samples=data$meta[[ data$sample_id ]],
      PTID=factor(data$meta[[ data$individual_id ]])
    )
    
    dt <- merge(dt, ptid_dt, all.x=TRUE, by="Samples")
    
    ## Sum the counts over PTID
    summed <- dt[, sum(Counts, na.rm=TRUE), by=PTID]
    output <- setNames(summed$V1, summed$PTID)
    return( output[ order(names(output)) ] )
    
  } else {
    return( setNames(dt$Counts, dt$Samples) )
  }
  
  
}
