\name{MakeComparable}
\alias{MakeComparable}
\title{Make a list of COMPASSResults comparable for display as heatamps}
\usage{
MakeComparable(X, subset, ...)
}
\arguments{
  \item{X}{An list of \code{COMPASSResult} objects}

  \item{subset}{An \R expression, evaluated within the
  metadata, used to determine which individuals should be
  kept.}

  \item{...}{Optional arguments passed to \code{pheatmap}.}
}
\description{
This function makes a list of COMPASSResults comparable for
display as heatmaps. It has some inefficiencies as it looks
at all pairwise comparisons.
}
\examples{
\dontrun{
#X is a list of COMPASSResult objects
X <- MakeComparable(X)
# These should show the same categories.
plot(X[[1]], remove_unexpressed_categories=FALSE)
plot(X[[2]], remove_unexpressed_categories=FALSE)
}
}

