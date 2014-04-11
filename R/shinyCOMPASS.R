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
##' @param facet1,facet2,facet3 Default values for facets in the Shiny app.
##'   Each should be the name of a single vector in the metadata.
##' @param main A title to give to the heatmap and subset histogram plots.
##' @param launch Boolean; if \code{TRUE} we launch the Shiny application.
##'   Otherwise, the user can launch it manually by navigating to the directory
##'   \code{dir} and running \code{shiny::runApp()}.
##' @param ... Optional arguments passed to \code{shiny::runApp}.
##' @seealso \code{\link{shinyCOMPASSDeps}}, for identifying packages that you
##'   need in order to run the Shiny application.
##' @export
##' @examples 
##' if (interactive()) {
##'   oldOpt <- getOption("example.ask")
##'   options(example.ask=FALSE)
##'   on.exit( options(example.ask=oldOpt) )
##'   shinyCOMPASS(CR)
##'   options(example.ask=TRUE)
##' }
shinyCOMPASS <- function(x, dir=NULL, meta.vars, facet1="None", facet2="None", 
  facet3="None", main="Heatmap of Ag-Specificity Posterior Probabilities", launch=TRUE, ...) {
  
  if (length(shinyCOMPASSDeps(verbose=FALSE))) {
    message("Error: Can't run the Shiny application as required packages ",
      "are missing. Instructions follow:\n\n")
    return(shinyCOMPASSDeps())
  }
  
  if (!require(shiny)) {
    stop("You must have 'shiny' installed to run the Shiny application -- try 'install.packages(\"shiny\")'.",
      call.=FALSE)
  }
  
  if (!inherits(x, "COMPASSResult"))
    stop("'shinyCOMPASS' can only be called on a COMPASSResult object", call.=FALSE)
  
  if (inherits(x, "SimpleCOMPASSResult")) {
    stop("'shinyCOMPASS' cannot be called on fits generated through 'SimpleCOMPASS'.",
      call.=FALSE)
  }
  
  message("Preparing data for the Shiny application, please wait a moment...")
  
  if (is.null(dir)) {
    dir <- file.path( tempdir(), "shinyCOMPASS" )
    on.exit(unlink(dir, recursive=TRUE))
  }
  
  dir <- path.expand(dir)
  
  ## Keep only the metadata variables specified
  iid <- x$data$individual_id
  sid <- x$data$sample_id
  
  if (!missing(meta.vars)) {
    
    x$data$meta <- x$data$meta[ 
      names(x$data$meta) %in% c(iid, sid, meta.vars)
    ]
    
    x$orig$meta <- x$orig$meta[
      names(x$orig$meta) %in% c(iid, sid, meta.vars)
    ]
  }
  
  ## Add the default facets
  x$facet1 <- facet1
  x$facet2 <- facet2
  x$facet3 <- facet3
  x$main   <- main
  
  ## Check the facets
  .check_facet <- function(facet) {
    if (facet != "None" && !(facet %in% names(x$orig$meta))) {
      stop("facet '", facet, "' is not available in the metadata")
    }
  }
  invisible(lapply(c(facet1, facet2, facet3), .check_facet))
  
  ## Copy the Shiny infrastructure files to the directory
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  file.copy(
    file.path( system.file(package="COMPASS"), "shiny/." ),
    file.path(dir),
    recursive=TRUE
  )
  
  ## Copy the data
  dir.create( file.path(dir, "data"), showWarnings=FALSE )
  saveRDS(x, file=file.path(dir, "data", "data.rds"))
  
  message("The files necessary for launching the COMPASS Shiny application have ", 
    "been copied to '", dir, "'.")
  
  if (launch) {
    message("Starting the Shiny application...")
    runApp(file.path(dir), ...)
  } else {
    file.path(dir)
  }
  
}
