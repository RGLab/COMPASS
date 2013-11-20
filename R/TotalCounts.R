##' Compute Total Cell Counts by Stimulation
##' 
##' @param data A \code{COMPASSContainer}.
##' @param subset An expression, evaluated within the metadata, defining
##'   the subset of \code{data} over which the counts are computed.
##' @param aggregate Boolean; if \code{TRUE} we sum over the individual,
##'   to get total counts across samples for each individual.
##' @export
TotalCounts <- function(data, subset, aggregate=TRUE) {
  
  if (!inherits(data, "COMPASSContainer"))
    stop("'data' must be an object of class 'COMPASSContainer'")
  
  if (!is.call(subset)) {
    subset_call <- match.call()$subset
  } else {
    subset_call <- subset
  }
  
  keep <- eval(subset_call, envir=data$meta)
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
      PTID=data$meta[[ data$individual_id ]]
    )
    
    dt <- merge(dt, ptid_dt, all.x=TRUE, by="Samples")
    
    ## Sum the counts over PTID
    output <- dt[, sum(Counts), by=PTID]
    return( setNames(output$V1, output$PTID) )
    
  } else {
    
    return (setNames(dt$Counts, dt$Samples))
    
  }
  
  
}
