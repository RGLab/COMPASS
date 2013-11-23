\name{FunctionalityScore}
\alias{FunctionalityScore}
\alias{FunctionalityScore.COMPASSResult}
\alias{FunctionalityScore.default}
\title{Compute the Functionality Score for each subject fit in a COMPASS model}
\usage{
  FunctionalityScore(x)

  \method{FunctionalityScore}{COMPASSResult} (x)

  \method{FunctionalityScore}{default} (x)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}, as
  returned by \code{\link{COMPASS}}. Alternatively, a
  matrix of functionality scores.}
}
\value{
  A numeric vector of functionality scores.
}
\description{
  Computes the functionality score for each observation
  from the gamma matrix of a COMPASS model fit.
}

