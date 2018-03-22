##' List Shiny Dependencies
##'
##' This function can be used to identify the packages still needed in order
##' to launch the Shiny app.
##'
##' @param verbose Boolean; if \code{TRUE} we print installation instructions
##'   to the screen.
##' @export
##' @examples
##' shinyCOMPASSDeps()
shinyCOMPASSDeps <- function(verbose=TRUE) {

  message <- function(...) {
    if (verbose) {
      return( base::message(...) )
    } else {
      return( invisible(NULL) )
    }
  }

  ## A flag that tracks whether an installation of one or more packages is
  ## necessary
  flag_install_required <- FALSE

  if (!requireNamespace("devtools", quietly = TRUE)) {
    flag_install_required <- TRUE
    message("You will need to install 'devtools' in order to download some of ",
      "the necessary packages: try install.packages('devtools').")
  }

  needed_pkgs <- list(
    CRAN=c("hexbin", "scales", "gridExtra", "ggplot2", "shiny", "reshape2",
      "data.table", "gtools", "stringr", "RColorBrewer"),
    KevinGitHub=c("data.table.extras"),
    WinstonGitHub="shinyGridster"
  )

  installed_pkgs <- rownames(installed.packages())
  CRAN_install <- setdiff(needed_pkgs$CRAN, installed_pkgs)

  if (length(CRAN_install)) {
    flag_install_required <- TRUE
    message("The following package(s) need to be installed from CRAN:\n\t",
      paste(CRAN_install, collapse=", "))
  }

  kevin_install <- setdiff(needed_pkgs$KevinGitHub, installed_pkgs)
  if (length(kevin_install)) {
    flag_install_required <- TRUE
    message("The following package(s) from GitHub need to be installed:\n\t",
      paste(kevin_install, collapse=", "), "\n\n",
      "These can be installed with 'devtools::install_github(\"kevinushey/<package_name>\")'"
    )
  }

  winston_install <- setdiff(needed_pkgs$WinstonGitHub, installed_pkgs)
  if (length(winston_install)) {
    flags_install_required <- TRUE
    message("The following pacakge(s) from GitHub need to be installed:\n\t",
      paste(winston_install, collapse=", "), "\n\n",
      "These can be installed with 'devtools::install_github(\"wch/<package_name>\")'"
    )
  }
  if (!flag_install_required) {
    message("Your system is ready to run the COMPASS Shiny application!\n",
      "Try calling 'shinyCOMPASS' on a COMPASS fit object to start the Shiny application.")
    return (invisible(NULL))
  }

  return( list(
    CRAN=CRAN_install,
    KevinGitHub=kevin_install,
    WinstonGitHub=winston_install
  ) )
}
