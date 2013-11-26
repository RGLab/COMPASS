\name{TotalCounts}
\alias{TotalCounts}
\title{Compute Total Cell Counts by Stimulation}
\usage{
  TotalCounts(data, subset, aggregate = TRUE)
}
\arguments{
  \item{data}{A \code{COMPASSContainer}.}

  \item{subset}{An expression, evaluated within the
  metadata, defining the subset of \code{data} over which
  the counts are computed.}

  \item{aggregate}{Boolean; if \code{TRUE} we sum over the
  individual, to get total counts across samples for each
  individual.}
}
\description{
  Compute Total Cell Counts by Stimulation
}

