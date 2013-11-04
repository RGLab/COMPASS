##' Start a Shiny Application for Visualizing COMPASS Results
##' 
##' This function takes a \code{COMPASSResult} object, and generates
##' a local Shiny application for visualizing the results.
##' 
##' @param x An object of class \code{COMPASSResult}.
##' @param dir A location to write out the \code{.rds} files that
##'   will be loaded and used by the Shiny application.
##' @param stimulated The name of the positive stimulation, as available
##'   in the metadata column specified by the \code{stimulation_id}.
##' @param unstimulated The name of the negative stimulation, as
##'   available in the metadata column specified by the \code{stimulation_id}.
##' @export
shinyCOMPASS <- function(x, dir=NULL,
  stimulated, unstimulated) {
  
  if (is.null(dir)) {
    dir <- tempdir()
    on.exit(unlink(dir, recursive=TRUE))
  }
  
  if (!inherits(x, "COMPASSResult"))
    stop("shinyCOMPASS can only be called on a COMPASSResult object")
  
  call <- match.call()
  
  ## We must pre-process the data to get it into Shiny
  dat <- x$orig$data
  meta <- x$orig$meta
  counts <- x$orig$counts
  Mgamma <- x$fit$mean_gamma
  categories <- x$fit$categories
  markers <- unname(colnames(x$orig$data[[1]]))
  
  colnames(Mgamma) <- apply(categories[, -ncol(categories)], 1, function(x) {
    paste0( swap(x, c("0", "1"), c("-", "+")), collapse="")
  })
  
  ## Compute the joint distribution of counts
  combos <- discrete_combinations( length(markers) )
  
  d_counts <- cell_counts(dat, combos)
  rownames(d_counts) <- names(dat)
  colnames(d_counts) <- sapply(combos, function(x) {
    paste0( swap(x, c(1:6, -1:-6), c(rep("+", 6), rep("-", 6))), collapse="" )
  })
  
  ## Only grab the metadata for which we have samples
  meta_ <- data.table(meta[ meta[[ x$data$sample_id ]] %in% names(dat), ])
  setnames(meta_, "name", "Sample")
  
  ## Melt and Merge
  m <- data.table(melt_(d_counts))
  setnames(m, c("Sample", "Cytokine", "Counts"))
  
  ## Generate columns for the markers
  ## 1 == on, 0 == off
  m[, (markers) := {
    transpose_list(lapply(strsplit( as.character(Cytokine), "", fixed=TRUE), function(x) {
      as.numeric(swap(x, c("+", "-"), c("1", "0")))
    }))
  }]
  
  ## Merge in the metadata
  setkey(m, Sample)
  setkey(meta_, Sample)
  d <- meta_[m]
  
  ## Compute the Log Fold Change
  stim <- x$orig$stimulation_id
  d[, LogFoldChange := log2(
    (Counts+1) / (Counts[ .SD[[stim]] == unstimulated ] + 1)
  ), by=list(PTID, Cytokine)]
  
  ## Merge in the MeanGamma
  ## we want to fill Mgamma with zeroes for the cytokines not included
  n_markers <- ncol(x$orig$data[[1]])
  Mgamma <- as.data.frame(Mgamma)
  combos <- do.call( paste0, 
    do.call( expand.grid, replicate(n_markers, c("-", "+"), simplify=FALSE) ) 
  )
  for (combo in combos) {
    if (!(combo %in% colnames(Mgamma))) {
      Mgamma[combo] <- 0
    }
  }
  
  Mgamma <- Mgamma[ order(colnames(Mgamma)) ]
  Mgamma <- as.matrix(Mgamma)
  M <- data.table(melt_(Mgamma))
  setnames(M, c("PTID", "Cytokine", "MeanGamma"))
  
  setkeyv(d, c("PTID", "Cytokine"))
  setkeyv(M, c("PTID", "Cytokine"))
  d <- M[d]
  
  dir.create(dir, showWarnings=FALSE)
  dir.create( file.path(dir, "data") )
  saveRDS(d, file=file.path(dir, "data", "data.rds"))
  
  ## Copy the ui.R, server.R used for starting the Shiny app
  file.copy( 
    system.file("shiny", "ui.R", package="COMPASS"),
    file.path(dir, "ui.R")
  )
  
  file.copy(
    system.file("shiny", "server.R", package="COMPASS"),
    file.path(dir, "server.R")
  )
    
  
}