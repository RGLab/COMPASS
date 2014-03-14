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
shinyCOMPASS <- function(x, dir=NULL, meta.vars, obfuscate=FALSE) {
  
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
  
  message("Starting the Shiny application...")
  runApp( file.path(dir) )
  
}
