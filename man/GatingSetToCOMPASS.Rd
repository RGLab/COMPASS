\name{GatingSetToCOMPASS}
\alias{GatingSetToCOMPASS}
\title{Extract A flowWorkspace Node to a COMPASSContainer}
\usage{
  GatingSetToCOMPASS(gs, node, meta, individual_id,
    sample_id)
}
\arguments{
  \item{gs}{A \code{GatingSet} object.}

  \item{node}{The name of the node.}

  \item{meta}{A \code{data.frame} of meta data. If it is
  missing, we attempt to pull the metadata from the
  \code{GatingSet}, and throw an error if it is
  unavailable.}

  \item{individual_id}{The name of the vector in
  \code{meta} that denotes the individuals from which
  samples were drawn.}

  \item{sample_id}{The name of the vector in \code{meta}
  that denotes the samples.  This vector should contain all
  of the names in the \code{data} input.}
}
\description{
  This function can be used to take gated data from a
  \code{flowWorkspace} object (e.g., a \code{GatingSet}),
  and retrieve it as a \code{COMPASSContainer}, suitable
  for ##' modelling with \code{COMPASS}.
}

