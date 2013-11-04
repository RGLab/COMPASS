\name{GatingSetToCOMPASS}
\alias{GatingSetToCOMPASS}
\title{Extract A flowWorkspace Node to a COMPASSContainer}
\usage{
  GatingSetToCOMPASS(gs, node, children, meta,
    individual_id, sample_id)
}
\arguments{
  \item{gs}{A \code{GatingSet} object.}

  \item{node}{The name of the node. We extract all children
  of that node, under the assumption that the markers of
  interest have been gated at the chosen node.}

  \item{children}{A character vector, or list, giving the
  names of the children of \code{node} that we wish to
  extract. If unspecified, we take all children of that
  node. If the character vector passed is named, we map
  from value to name; e.g. if a node is named \code{154}
  but we prefer the name \code{CD154+}, we could pass
  \code{`CD154+`="154"}.}

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
\seealso{
  \code{\link{COMPASSContainer}}
}

