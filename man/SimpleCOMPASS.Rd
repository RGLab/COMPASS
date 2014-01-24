\name{SimpleCOMPASS}
\alias{SimpleCOMPASS}
\title{Fit the discrete COMPASS Model}
\usage{
SimpleCOMPASS(n_s, n_u, meta, individual_id, sample_id, iterations = 10000L,
  replications = 8, verbose = TRUE)
}
\arguments{
  \item{n_s}{The cell counts for stimulated cells.}

  \item{n_u}{The cell counts for unstimulated cells.}

  \item{meta}{A \code{data.frame} of metadata, describing
  the individuals in the experiment. Each row in
  \code{meta} should correspond to a row in \code{data}.
  There should be one row for each sample; i.e., one row
  for each element of \code{data}.}

  \item{individual_id}{The name of the vector in
  \code{meta} that denotes the individuals from which
  samples were drawn.}

  \item{sample_id}{The name of the vector in \code{meta}
  that denotes the samples.  This vector should contain all
  of the names in the \code{data} input.}

  \item{iterations}{The number of iterations (per
  'replication') to perform.}

  \item{replications}{The number of 'replications' to
  perform. In order to conserve memory, we only keep the
  model estimates from the last replication.}

  \item{verbose}{Boolean; if \code{TRUE} we output progress
  information.}
}
\value{
A \code{list} with class \code{COMPASSResult} with two
components, the \code{fit} containing parameter estimates
and parameter acceptance rates, and \code{data} containing
the generated data used as input for the model.
}
\description{
This function fits the \code{COMPASS} model from a
user-provided set of stimulated / unstimulated matrices.
}

