## debug
if (FALSE) {
  library(COMPASS)
  x <- readRDS("data/RV144_CD4_results_discrete.rds")
  dir <- "shinyTest"
  stimulated <- "92TH023 Env"
  unstimulated <- "negctrl 1"
  shinyCOMPASS(x, dir, stimulated, unstimulated)
}

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
  
  if (missing(stimulated)) {
    stop("'stimulated' must be supplied; it is necessary for calculations of ",
      "Log Fold Change and other variables used in the Shiny application")
  }
  
  if (missing(stimulated)) {
    stop("'unstimulated' must be supplied; it is necessary for calculations of ",
      "Log Fold Change and other variables used in the Shiny application")
  }
  
  message("Preparing data for the Shiny application, please wait a moment...")
  
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
  iid <- x$orig$individual_id
  sid <- x$orig$sample_id
  stid <- x$orig$stimulation_id
  
  colnames(Mgamma) <- apply(categories[, -ncol(categories)], 1, function(x) {
    paste0( COMPASS:::swap(x, c("0", "1"), c("-", "+")), collapse="")
  })
  
  ## Compute the joint distribution of counts
  combos <- COMPASS:::discrete_combinations( length(markers) )
  
  d_counts <- COMPASS:::cell_counts(dat, combos)
  rownames(d_counts) <- names(dat)
  colnames(d_counts) <- sapply(combos, function(x) {
    paste0( swap(x, c(1:6, -1:-6), c(rep("+", 6), rep("-", 6))), collapse="" )
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
    (Counts+1) / (Counts[ .SD[[stim]] == unstimulated ] + 1)
  ), by=c(iid, "Marker")]
  
  ## Compute the difference in proportions
  d[, PropDiff := Proportion - Proportion[ .SD[[stim]] == unstimulated ],
    by=c(iid, "Marker")]
  
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
  
  dir.create(dir, showWarnings=FALSE)
  
  ## Copy the Shiny infrastructure to the directory
  file.copy(
    file.path(system.file(package="COMPASS"), "shiny"),
    file.path(dir),
    recursive=TRUE
  )
  
  ## Copy the data
  dir.create( file.path(dir, "shiny", "data"), showWarnings=FALSE )
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
  saveRDS(output, file=file.path(dir, "shiny", "data", "data.rds"))
  
  message("Starting the Shiny application...")
  shiny::runApp( file.path(dir, "shiny") )
  
}
