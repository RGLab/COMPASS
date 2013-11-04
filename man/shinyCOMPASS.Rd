\name{shinyCOMPASS}
\alias{shinyCOMPASS}
\title{Start a Shiny Application for Visualizing COMPASS Results}
\usage{
  shinyCOMPASS(x, dir = NULL, stimulated, unstimulated)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}.}

  \item{dir}{A location to write out the \code{.rds} files
  that will be loaded and used by the Shiny application.}

  \item{stimulated}{The name of the positive stimulation,
  as available in the metadata column specified by the
  \code{stimulation_id}.}

  \item{unstimulated}{The name of the negative stimulation,
  as available in the metadata column specified by the
  \code{stimulation_id}.}
}
\description{
  This function takes a \code{COMPASSResult} object, and
  generates a local Shiny application for visualizing the
  results.
}

