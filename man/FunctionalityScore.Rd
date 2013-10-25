\name{FunctionalityScore}
\alias{FunctionalityScore}
\title{Compute the Functionality Score for each subject fit in a COMPASS model}
\usage{
  FunctionalityScore(CR)
}
\arguments{
  \item{CR}{An object of class \code{COMPASSResult}, as
  returned by \code{\link{COMPASS}}.}
}
\value{
  A numeric vector of functionality scores.
}
\description{
  Computes the functionality score for each observation
  from the gamma matrix of a COMPASS model fit.
}

