\name{TotalCellCounts}
\alias{TotalCellCounts}
\title{Compute Total Cell Counts}
\usage{
TotalCellCounts(data, subset, aggregate = TRUE)
}
\arguments{
  \item{data}{A \code{COMPASSContainer}.}

  \item{subset}{An expression, evaluated within the
  metadata, defining the subset of \code{data} over which
  the counts are computed. If left unspecified, the counts
  are computed over all samples.}

  \item{aggregate}{Boolean; if \code{TRUE} we sum over the
  individual, to get total counts across samples for each
  individual.}
}
\description{
This function is used to compute total cell counts, per
individual, from a \code{COMPASSContainer}.
}
\examples{
TotalCellCounts(CC, trt == "Treatment")
TotalCellCounts(CC, trt == "Control")
TotalCellCounts(CC)
}

