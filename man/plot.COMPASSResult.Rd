\name{plot.COMPASSResult}
\alias{plot}
\alias{plot.COMPASSResult}
\title{Plot a COMPASSResult}
\usage{
\method{plot}{COMPASSResult}(x, y, subset,
  remove_unexpressed_categories = TRUE, minimum_dof = 1,
  maximum_dof = Inf, row_annotation,
  palette = colorRampPalette(brewer.pal(10, "Purples"))(20),
  show_rownames = FALSE, show_colnames = FALSE, measure = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}.}

  \item{y}{This argument gets passed to
  \code{row_annotation}, if \code{row_annotation} is
  missing. It can be used to group rows (individuals) by
  different conditions as defined in the metadata.}

  \item{subset}{An \R expression, evaluated within the
  metadata, used to determine which individuals should be
  kept.}

  \item{remove_unexpressed_categories}{Boolean, if
  \code{TRUE} we remove any unexpressed categories.}

  \item{minimum_dof}{The minimum degree of functionality
  for the categories to be plotted.}

  \item{maximum_dof}{The maximum degree of functionality
  for the categories to be plotted.}

  \item{row_annotation}{A vector of names, pulled from the
  metadata, to be used for row annotation.}

  \item{palette}{The colour palette to be used.}

  \item{show_rownames}{Boolean; if \code{TRUE} we display
  row names (ie, the individual ids).}

  \item{show_colnames}{Boolean; if \code{TRUE} we display
  column names (ie, the column name associated with a
  cytokine; typically not needed)}

  \item{measure}{Optional. By default, we produce a heatmap
  of the mean gammas produced in a model fit. We can
  override this by supplying a matrix of suitable dimension
  as well; these can be generated with the
  \code{Posterior*} functions -- see
  \code{\link{Posterior}} for examples.}

  \item{...}{Optional arguments passed to \code{pheatmap}.}
}
\description{
This function can be used to visualize the mean probability
of response -- that is, the probability that there is a
difference in response between samples subjected to the
'treatment' condition, and samples subjected to the
'control' condition.
}
\examples{
## visualize the mean probability of reponse
plot(CR)

## visualize the proportion of cells belonging to a category
plot(CR, measure=PosteriorPs(CR))
}

