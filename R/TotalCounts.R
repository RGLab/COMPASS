##' Compute Total Cell Counts by Stimulation
##' 
##' @param data A \code{COMPASSContainer}.
##' @param stimulation A stimulation, or set of stimulations, expressed as a
##'   boolean combination.
##' @export
TotalCounts <- function(data, stimulation) {
  
  if (!inherits(data, "COMPASSContainer"))
    stop("'data' must be an object of class 'COMPASSContainer'")
  
  keep <- data$meta[[ data$stimulation_id ]] %like% stimulation
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
  
  ## Merge in the PTID information
  ptid_dt <- data.table(
    Samples=data$meta[[ data$sample_id ]],
    PTID=data$meta[[ data$individual_id ]]
  )
  
  dt <- merge(dt, ptid_dt, all.x=TRUE, by="Samples")
  
  ## Sum the counts over PTID
  output <- dt[, sum(Counts), by=PTID]
  return( setNames(output$V1, output$PTID) )
  
}
