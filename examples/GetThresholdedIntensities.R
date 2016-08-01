if (require("flowWorkspace")) {

  ## Generate an example GatingSet that could be used with COMPASS
  ## We then pull out the 'data' and 'counts' components that could
  ## be used within a COMPASSContainer

  n <- 10 ## number of samples
  k <- 4 ## number of markers

  sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
  iid_vec <- rep_len( paste0("iid_", 1:(n/10) ), n ) ## individual ids
  marker_names <- c("TNFa", "IL2", "IL4", "IL6")

  ## Generate n sets of 'flow' data -- a list of matrices, each row
  ## is a cell, each column is fluorescence intensities on a particular
  ## channel / marker
  data <- replicate(n, {
    nrow <- round(runif(1) * 1E4 + 1000)
    ncol <- k
    vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
    output <- matrix(vals, nrow, ncol)
    colnames(output) <- marker_names
    return(output)
  })
  names(data) <- sid_vec

  ## Put it into a GatingSet
  fs <- flowSet( lapply(data, flowFrame) )
  gs <- GatingSet(fs)

  ## Add some dummy metadata
  meta <- pData(gs)
  meta$PTID <- 1:10
  pData(gs) <- meta

  gate <- rectangleGate( list(TNFa=c(-Inf,Inf)))
  add(gs, gate, parent="root", name="dummy")

  ## Add dummy gate

  ## Make some gates, and apply them
  invisible(lapply(marker_names, function(marker) {
    .gate <- setNames( list( c( rexp(1, runif(1, 1E-5, 1E-3)), Inf) ), marker )
    gate <- rectangleGate(.gate=.gate)
    add(gs, gate, parent="dummy", name=paste0(marker, "+"))
  }))

  recompute(gs)

  ## Map node names to channel names
  map=list(
    "TNFa+"="TNFa",
    "IL2+"="IL2",
    "IL4+"="IL4",
    "IL6+"="IL6"
  )

  ## Pull out the data as a COMPASS-friendly dataset
  node <- "dummy"
  map <- map
  system.time(
    output <- GetThresholdedIntensities(gs, "dummy", map)
  )

  system.time(
    output <- COMPASSContainerFromGatingSet(gs, "dummy", individual_id="PTID")
  )

  str(output)

}
