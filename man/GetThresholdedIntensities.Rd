\name{GetThresholdedIntensities}
\alias{GetThresholdedIntensities}
\title{Extract Thresholded Intensities from a GatingSet}
\usage{
GetThresholdedIntensities(gs, node, map)
}
\arguments{
  \item{gs}{A \code{GatingSet} or \code{GatingSetList}.}

  \item{node}{The name, or path, of a single node in a
  \code{GatingSet} / \code{GatingSetList}.}

  \item{map}{A \code{list}, mapping node names to markers.}
}
\value{
A \code{list} with two components:

\item{data}{A \code{list} of thresholded intensity
measures.} \item{counts}{A named vector of total cell
counts at the node \code{node}.}
}
\description{
This function extracts thresholded intensities for children
of a node \code{node}, as specified through the \code{map}
argument.
}
\details{
\code{map} should be an \R \code{list}, mapping node names
(as specified in the gating hierarchy of the gating set) to
channel names (as specified in either the \code{desc} or
\code{name} columns of the parameters of the associated
\code{flowFrame}s in the \code{GatingSet}).
}
\examples{
\dontrun{
  GetThresholdedIntensities(gs, "4+", list("INFg+"="IFNg"))
}
}

