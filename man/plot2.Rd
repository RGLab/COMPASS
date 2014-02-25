\name{plot2}
\alias{plot2}
\title{Plot a pair of COMPASSResults}
\usage{
plot2(x, y, row_annotation = NULL, remove_unexpressed_categories = TRUE,
  minimum_dof = 1, maximum_dof = Inf, subset, palette = NA,
  show_rownames = FALSE, show_colnames = FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}.}

  \item{y}{An object of class \code{COMPASSResult}.}

  \item{row_annotation}{A vector of names, pulled from the
  metadata, to be used for row annotation.}

  \item{remove_unexpressed_categories}{Boolean, if
  \code{TRUE} we remove any unexpressed categories.}

  \item{minimum_dof}{The minimum degree of functionality
  for the categories to be plotted.}

  \item{maximum_dof}{The maximum degree of functionality
  for the categories to be plotted.}

  \item{subset}{An \R expression, evaluated within the
  metadata, used to determine which individuals should be
  kept.}

  \item{palette}{The colour palette to be used.}

  \item{show_rownames}{Boolean; if \code{TRUE} we display
  row names (ie, the individual ids).}

  \item{show_colnames}{Boolean; if \code{TRUE} we display
  column names (ie, the column name associated with a
  cytokine; typically not needed)}

  \item{...}{Optional arguments passed to \code{pheatmap}.}
}
\value{
The plot as a \code{grid} object (\code{grob}). It can be
redrawn with e.g. \code{grid::grid.draw()}.
}
\description{
This function can be used to visualize the mean probability
of response -- that is, the probability that there is a
difference in response between samples subjected to the
'treatment' condition, and samples subjected to the
'control' condition.
}

