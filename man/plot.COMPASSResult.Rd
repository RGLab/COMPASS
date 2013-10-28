\name{plot.COMPASSResult}
\alias{plot.COMPASSResult}
\title{Plot a Heatmap of the Mean Gammas}
\usage{
  \method{plot}{COMPASSResult} (x, y, subset,
    row_annotation, ...)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}.}

  \item{y}{This argument gets passed to
  \code{row_annotation}, if \code{row_annotation} is
  missing.}

  \item{subset}{An \R expression, evaluated within the
  metadata, used to determine which individuals should be
  kept.}

  \item{row_annotation}{A vector of names, pulled from the
  metadata, to be used for row annotation.}

  \item{...}{Optional arguments passed to \code{pheatmap}.}
}
\description{
  This function plots the mean gammas.
}

