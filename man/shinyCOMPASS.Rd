\name{shinyCOMPASS}
\alias{shinyCOMPASS}
\title{Start a Shiny Application for Visualizing COMPASS Results}
\usage{
shinyCOMPASS(x, dir = NULL, meta.vars, obfuscate = FALSE)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}.}

  \item{dir}{A location to write out the \code{.rds} files
  that will be loaded and used by the Shiny application.}

  \item{meta.vars}{A character vector of column names that
  should be used for potential facetting in the Shiny app.
  By default, we take all metadata variables; you may want
  to limit this if you know certain variables are not of
  interest.}

  \item{obfuscate}{Boolean; if \code{TRUE} we replace the
  patient IDs used with a randomly generated set of IDs.
  This is useful if you wish to display the data you are
  using externally but don't wish to make available the IDs
  used.}
}
\description{
This function takes a \code{COMPASSResult} object, and
generates a local Shiny application for visualizing the
results.
}
\examples{
if (interactive()) {
  oldOpt <- getOption("example.ask")
  options(example.ask=FALSE)
  on.exit( options(example.ask=oldOpt) )
  shinyCOMPASS(CR)
  options(example.ask=TRUE)
}
}
\seealso{
\code{\link{shinyCOMPASSDeps}}, for identifying packages
that you need in order to run the Shiny application.
}

