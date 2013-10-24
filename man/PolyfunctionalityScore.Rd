\name{PolyfunctionalityScore}
\alias{PolyfunctionalityScore}
\title{Compute the Polyfunctionality Score for each subject fit in a COMPASS model}
\usage{
  PolyfunctionalityScore(fit,
    normalization = c("all", "observed"))
}
\arguments{
  \item{fit}{is a \code{COMPASSFit} object from either the
  discrete or continuous COMPASS model.}

  \item{normalization}{a \code{character} vector specifying
  how the score is to be normalized. Either using
  \code{"all"} possible categories, or just the
  \code{"observed"} categories. Defaults to \code{"all"}
  categories.}
}
\value{
  a \code{vector} of polyfunctionality scores.
}
\description{
  Computes the Polyfunctionality score for each observation
  from the gamma matrix of a COMPASS model fit.
}

