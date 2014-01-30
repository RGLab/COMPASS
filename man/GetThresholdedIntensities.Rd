\name{GetThresholdedIntensities}
\alias{GetThresholdedIntensities}
\title{Extract Thresholded Counts from a GatingSet}
\usage{
GetThresholdedIntensities(gs, node, map)
}
\arguments{
  \item{gs}{A \code{GatingSet} or \code{GatingSetList}.}

  \item{node}{The name, or path, of a single node in a
  \code{GatingSet} / \code{GatingSetList}.}

  \item{map}{A \code{list}, mapping node names to markers.}
}
\description{
This function extracts thresholded intensities for children
of a node \code{node}, as specified through the \code{map}
argument.
}
\details{
\code{map} should be a list mapping node names, as
specified in the gating hierarchy of the gating set, to
channel names, as specified in either the \code{desc} or
\code{name} columns of the parameters of the associated
\code{flowFrame}s in the \code{GatingSet}.
}
\examples{
\dontrun{
  GetThresholdedIntensities(gs, "4+")
}
}

