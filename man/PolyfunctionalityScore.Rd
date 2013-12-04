\name{PolyfunctionalityScore}
\alias{PolyfunctionalityScore}
\alias{PolyfunctionalityScore.COMPASSResult}
\alias{PolyfunctionalityScore.default}
\title{Compute the Polyfunctionality Score for each subject fit in a COMPASS model}
\usage{
  PolyfunctionalityScore(x, degree, normalization)

  \method{PolyfunctionalityScore}{COMPASSResult} (x,
    degree, normalization = "all")

  \method{PolyfunctionalityScore}{default} (x, degree,
    normalization)
}
\arguments{
  \item{x}{An object of class \code{COMPASSResult}, as
  returned by \code{\link{COMPASS}}.}

  \item{degree}{TODO}

  \item{normalization}{A \code{character} vector specifying
  how the score is to be normalized. Either using
  \code{"all"} possible categories, or just the
  \code{"observed"} categories. Defaults to \code{"all"}
  categories.}
}
\value{
  A numeric vector of polyfunctionality scores.
}
\description{
  Computes the Polyfunctionality score for each observation
  from the gamma matrix of a COMPASS model fit.
}

