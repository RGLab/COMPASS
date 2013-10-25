\name{GatingHierarchyToCOMPASS}
\alias{GatingHierarchyToCOMPASS}
\title{Extract A flowWorkspace Node to a COMPASSContainer}
\usage{
  GatingHierarchyToCOMPASS(gh, node)
}
\arguments{
  \item{gh}{A \code{GatingHierarchy} object.}

  \item{node}{The name of the node.}
}
\description{
  This function can be used to take gated data from a
  \code{flowWorkspace} object (e.g., a
  \code{GatingHierarchy}), and retrieve it as a
  \code{COMPASSContainer}, suitable for ##' modelling with
  \code{COMPASS}.
}

