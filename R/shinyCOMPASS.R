##' Start a Shiny Application for Visualizing COMPASS Results
##' 
##' This function takes a \code{COMPASSResult} object, and generates
##' a local Shiny application for visualizing the results.
##' 
##' @param x An object of class \code{COMPASSResult}.
##' @param dir A location to write out the \code{.rds} files that
##'   will be loaded and used by the Shiny application.
##' @param meta.vars A character vector of column names that should be used
##'   for potential facetting in the Shiny app. By default, we take all
##'   metadata variables; you may want to limit this if you know certain
##'   variables are not of interest.
##' @param obfuscate Boolean; if \code{TRUE} we replace the patient IDs
##'   used with a randomly generated set of IDs. This is useful if you wish
##'   to display the data you are using externally but don't wish to make
##'   available the IDs used.
##' @export
shinyCOMPASS <- function(x, dir=NULL, meta.vars, obfuscate=FALSE) {
  
  if (!require(shiny)) {
    stop("You must have 'shiny' installed to run the Shiny application -- try 'install.packages(\"shiny\")'.",
      call.=FALSE)
  }
  
  if (!inherits(x, "COMPASSResult"))
    stop("'shinyCOMPASS' can only be called on a COMPASSResult object", call.=FALSE)
  
  message("Preparing data for the Shiny application, please wait a moment...")
  
  ## These are just R CMD check silencers
  ## They're all used as naked symbols in data.table code later
  Marker <- Proportion <- Counts <- LogFoldChange <- PropDiff <- Degree <-
    DOF <- ModelDiff <- ModelLogDiff <- NULL
  
  if (is.null(dir)) {
    dir <- tempdir()
    on.exit(unlink(dir, recursive=TRUE))
  }
  
  ## We must pre-process the data to get it into Shiny
  dat <- x$orig$data
  meta <- x$orig$meta
  counts <- x$orig$counts
  Mgamma <- x$fit$mean_gamma
  categories <- x$fit$categories
  markers <- unname(colnames(x$orig$data[[1]]))
  iid <- x$orig$individual_id
  sid <- x$orig$sample_id
  
  ## Obfuscate the ptids
  if (obfuscate) {
    
    n <- length( unique( meta[[iid]] ) )
    logn <- floor( log10(n) )
    .swap <- function(x) swap(x, from, to)
    
    from <- unique( meta[[ iid ]] )
    to <- paste("Individual", sprintf( paste0("%0", logn, "i"), 1:n))
    
    meta[[ iid ]] <- .swap( meta[[iid]] )
    rownames(Mgamma) <- .swap( rownames(Mgamma) )
    
  }
  
  ## Get the stimulated, unstimulated calls from the fit call
  call <- x$fit$call
  stimulated <- call[["treatment"]]
  unstimulated <- call[["control"]]
  
  ## Get the stimulation variable name from the call
  stid <- as.character(call[["treatment"]][[2]])
  if (!missing(meta.vars)) {
    meta <- meta[c(meta.vars, iid, sid, stid)]
  }
  
  colnames(Mgamma) <- apply(categories[, -ncol(categories)], 1, function(x) {
    paste0( swap(x, c("0", "1"), c("-", "+")), collapse="")
  })
  
  ## Compute the joint distribution of counts
  k <- length(markers)
  combos <- Combinations(k)
  
  d_counts <- CellCounts(dat, combos)
  rownames(d_counts) <- names(dat)
  colnames(d_counts) <- sapply(combos, function(x) {
    paste0( swap(x, c(1:k, -1:-k), c(rep("+", k), rep("-", k))), collapse="" )
  })
  
  ## Only grab the metadata for which we have samples
  meta_ <- data.table(meta[ meta[[ sid ]] %in% names(dat), ])
  
  ## Melt and Merge
  m <- data.table(melt_(d_counts))
  setnames(m, c(sid, "Marker", "Counts"))
  
  ## Generate columns for the markers
  ## 1 == on, 0 == off
  m[, (markers) := {
    transpose_list(lapply(strsplit( as.character(Marker), "", fixed=TRUE), function(x) {
      as.numeric(swap(x, c("+", "-"), c("1", "0")))
    }))
  }]
  
  ## Merge in the metadata
  setkeyv(m, sid)
  setkeyv(meta_, sid)
  d <- meta_[m]
  
  ## Compute the proportions
  
  ## First, merge in the total cell counts
  counts_dt <- data.table(
    names(counts),
    unname(counts)
  )
  
  setnames(counts_dt, c(sid, "TotalCellCounts"))
  setkeyv(counts_dt, sid)
  setkeyv(d, sid)
  d <- counts_dt[d]
  
  ## Next, compute proportions as Counts / TotalCellCounts
  d[, Proportion := Counts / TotalCellCounts]
  
  ## Compute the Log Fold Change
  stim <- x$orig$stimulation_id
  d[, LogFoldChange := log2(
    (Counts+1) / (mean(Counts[ eval(unstimulated) ]) + 1)
  ), by=c(iid, "Marker")]
  
  ## Compute the difference in proportions
  d[, PropDiff := Proportion - mean(Proportion[ eval(unstimulated) ]),
    by=c(iid, "Marker")]
  
  ## Compute the degree (== number of positive markers)
  d[, Degree := apply(.SD, 1, sum), .SDcols=markers]
  
  ## Compute the degree of functionality as the sum of 
  ## ps-pu over subsets of the same degree.
  d[, DOF := sum(PropDiff), by=c(iid, "Degree")]
  
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
  setnames(M, c(iid, "Marker", "MeanGamma"))
  
  setkeyv(d, c(iid, "Marker"))
  setkeyv(M, c(iid, "Marker"))
  d <- M[d]
  
  ## Merge the functionality, polyfunctionality scores
  scores <- data.table(
    FunctionalityScore=FunctionalityScore(x),
    PolyfunctionalityScore=PolyfunctionalityScore(x),
    PTID=names(FunctionalityScore(x))
  )
  setnames(scores, "PTID", iid)
  setkeyv(scores, iid)
  d <- scores[d]
  
  ## Merge in the model-estimated ps - pu, log(ps) - log(pu)
  cats <- apply(categories[-nrow(categories), -ncol(categories)], 1, function(x) {
    paste0( swap(x, c("0", "1"), c("-", "+")), collapse="")
  })
  
  model_ests <- data.table(
    ModelDiff=unlist(lapply(x$fit$posterior, "[[", "diff")),
    ModelLogDiff=unlist(lapply(x$fit$posterior, "[[", "logd")),
    PTID=rep( names(x$fit$posterior), each=length(x$fit$posterior[[1]][[1]]) ),
    category=cats
  )
  setnames(model_ests, c("PTID", "category"), c(iid, "Marker"))
  
  setkeyv(model_ests, c(iid, "Marker"))
  setkeyv(d, c(iid, "Marker"))
  d <- model_ests[d]
  d[is.na(ModelDiff), ModelDiff := 0]
  d[is.na(ModelLogDiff), ModelLogDiff := 0]
  
  ## Rename the Marker column
  d[, Marker := 
      unlist( lapply(strsplit(as.character(Marker), "", fixed=TRUE), function(x) {
        paste(markers, x, sep="", collapse="")
      }) )
    ]
  
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  file.copy(
    file.path( system.file(package="COMPASS"), "shiny/." ),
    file.path(dir),
    recursive=TRUE
  )
  
  ## Copy the data
  dir.create( file.path(dir, "data"), showWarnings=FALSE )
  output <- list(
    data=d, 
    meta=meta_, 
    cell_data=x$orig$data, 
    counts=x$orig$counts,
    markers=colnames(x$data$categories)[-ncol(x$data$categories)],
    COMPASS=x,
    stimulated=stimulated,
    unstimulated=unstimulated
  )
  saveRDS(output, file=file.path(dir, "data", "data.rds"))
  
  message("Starting the Shiny application...")
  runApp( file.path(dir) )
  
}
