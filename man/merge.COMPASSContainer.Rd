\name{merge.COMPASSContainer}
\alias{merge.COMPASSContainer}
\title{Merge Two COMPASSContainers}
\usage{
\method{merge}{COMPASSContainer}(x, y, ...)
}
\arguments{
  \item{x}{A \code{COMPASSContainer}.}

  \item{y}{A \code{COMPASSContainer}.}

  \item{...}{Optional arguments, currently ignored.}
}
\description{
This function merges two \code{COMPASSContainers}.
}
\examples{
## Chop the example COMPASSContainer into two, then merge it back together
CC1 <- subset(CC, trt == "Control")
CC2 <- subset(CC, trt == "Treatment")
merged <- merge(CC1, CC2)
res <- identical(CC, merge(CC1, CC2)) ## should return TRUE in this case
stopifnot( isTRUE(res) )
}

