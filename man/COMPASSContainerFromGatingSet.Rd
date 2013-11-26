\name{COMPASSContainerFromGatingSet}
\alias{COMPASSContainerFromGatingSet}
\title{Create a COMPASS Container from a GatingSet}
\usage{
  COMPASSContainerFromGatingSet(gs = NULL, node = NULL,
    filter.fun = NULL, individual_id = "PTID",
    sample_id = "name", stimulation_id = "Stim", mp = NULL,
    countFilterThreshold = 5000)
}
\arguments{
  \item{gs}{a \code{GatingSet} or \code{GatingSetList}}

  \item{node}{a \code{regular expression} to match a single
  node in the gating tree. If more than one node is
  matched, an error is thrown.}

  \item{filter.fun}{a \code{function} that does string
  substitution to clean up node names, i.e. turns a "CD4+"
  into a "CD4" to try and match against the
  \code{parameters} slot of the \code{flowFrames} in
  \code{gs}}

  \item{individual_id}{a \code{character} identifying the
  subject id column in the \code{gs} metadata}

  \item{sample_id}{a \code{character} idetifying the sample
  id column in the \code{gs} metadata.}

  \item{stimulation_id}{a \code{character} identifying the
  stimulation or treatment columnin the \code{gs}
  metadata.}

  \item{mp}{a \code{list} mapping node names to markers.
  This function tries to guess, but may fail. The user can
  override the guesswork.}
}
\description{
  This code expects a \code{GatingSet} or
  \code{GatingSetList}. It expects a regular expression for
  the node name (i.e. "/4\\+$" would match "/4+" in a node
  name with the plus sign at the end of the string. The
  user must supply the individual_id, sample_id, and
  stimulation_id, but they have default values suitable for
  the data we commonly see. Sometimes the child node names
  don't match the marker names exactly. This function will
  try to make some guesses about how to match these up. The
  \code{filter.fun} parameter is a function that does some
  regular expression string substitution to try and clean
  up the node names by removing various symobls that are
  often added to gates, \{+/-\}. The user can provide their
  own function to do string cleanup. Counts are extracted
  as well as metadata and single cell data, and these are
  fed into the COMPASSContainer constructor.
}
\details{
  There is likely not sufficient error checking.
}
\examples{
\dontrun{
  COMPASSContainerFromGatingSet(gatingset,"/4\\\\+$")
}
}
\seealso{
  \code{\link{COMPASSContainer}}
}

