\name{subset.COMPASSContainer}
\alias{subset.COMPASSContainer}
\title{Subset a COMPASSContainer}
\usage{
\method{subset}{COMPASSContainer}(x, subset, ...)
}
\arguments{
  \item{x}{A \code{COMPASSContainer}.}

  \item{subset}{A logical expression, evaluated within the
  metadata, indicating which entries to keep.}

  \item{...}{Optional arguents; currently unused.}
}
\description{
Use this function to subset a \code{COMPASSContainer}.
}
\examples{
subset(CC, iid == "iid_1")
}

