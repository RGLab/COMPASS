\docType{package}
\name{COMPASS}
\alias{COMPASS}
\alias{COMPASS-package}
\title{Fit the COMPASS Model}
\usage{
  COMPASS(data, treatment, control, subset = NULL,
    category_filter = function(x) colSums(x > 5) > 2,
    filter_lowest_frequency = 0,
    filter_specific_markers = NULL,
    model = c("discrete", "continuous"),
    iterations = 40000, replications = 8, verbose = TRUE,
    ...)
}
\arguments{
  \item{data}{An object of class \code{COMPASSData}.}

  \item{treatment}{An \R expression, evaluated within the
  metadata, that returns \code{TRUE} for those samples that
  should belong to the treatment group. For example, if the
  samples that received a positive stimulation were named
  \code{"92TH023 Env"} within a variable in \code{meta}
  called \code{Stim}, you could write \code{Stim ==
  "92TH023 Env"}.}

  \item{control}{An \R expression, evaluated within the
  metadata, that returns \code{TRUE} for those samples that
  should belong to the control group. See above for
  details.}

  \item{subset}{An expression used to subset the data. We
  keep only the samples for which the expression evaluates
  to \code{TRUE} in the metadata.}

  \item{category_filter}{A filter for the categories that
  are generated. This is a function that will be applied to
  the \emph{treatment counts} matrix generated from the
  intensities. Only categories meeting the
  \code{category_filter} criteria will be kept.}

  \item{filter_lowest_frequency}{A number specifying how
  many of the least expressed markers should be removed.}

  \item{filter_specific_markers}{Similar to
  \code{filter_lowest_frequency}, but lets you explicitly
  exclude markers.}

  \item{model}{A string denoting which model to fit;
  \code{"discrete"} indicates we should fit the discrete
  \code{COMPASS} model, while \code{"continuous"} indicates
  we should fit the continuous \code{COMPASS} model.}

  \item{iterations}{The number of iterations (per
  'replication') to perform.}

  \item{replications}{The number of 'replications' to
  perform. In order to conserve memory, we only keep the
  model estimates from the last replication.}

  \item{verbose}{Boolean; if \code{TRUE} we output progress
  information.}

  \item{...}{Other arguments; currently unused.}
}
\description{
  This function fits the \code{COMPASS} model.

  This package implements a model for the analysis of
  polyfunctionality in single-cell experiments.
}
\section{Category Filter}{
  The category filter is used to exclude categories
  (combinations of markers expressed for a particular cell)
  that are expressed very rarely. It is applied to the
  \code{treatment} \emph{counts} matrix, which is a
  \code{N} samples by \code{K} categories matrix. Those
  categories which are mostly unexpressed can be excluded
  here. For example, the default criteria,

  \code{category_filter=function(x) colSums(x > 5) > 2}

  indicates that we should only retain categories for which
  at least two samples had at least 5 cells expressing that
  particular combination of markers.
}
\seealso{
  \code{\link{COMPASSContainer}}, for constructing the data
  object required by \code{COMPASS}.
}

