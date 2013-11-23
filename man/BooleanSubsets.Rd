\name{BooleanSubsets}
\alias{BooleanSubsets}
\alias{BooleanSubsets.COMPASSContainer}
\alias{BooleanSubsets.default}
\title{Generate Boolean Subsets}
\usage{
  BooleanSubsets(x, as.matrix)

  \method{BooleanSubsets}{COMPASSContainer} (x,
    as.matrix = FALSE)

  \method{BooleanSubsets}{default} (x, as.matrix = FALSE)
}
\arguments{
  \item{x}{Either a \code{COMPASSContainer}, or a list of
  matrices.}

  \item{as.matrix}{Boolean; if \code{TRUE} we return
  results as a matrix; otherwise, we return the results as
  a list.}
}
\description{
  Generate all possible boolean subsets (that is, unique
  rows by index).
}
\details{
  \itemize{
  \item{\code{\link{BooleanSubsets.COMPASSContainer}}}{ for
  a \code{COMPASSContainer}},
  \item{\code{\link{BooleanSubsets.default}}}{ for the
  default method} }
}

