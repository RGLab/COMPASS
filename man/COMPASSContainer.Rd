\name{COMPASSContainer}
\alias{COMPASSContainer}
\title{Generate the Data Object used by COMPASS}
\usage{
  COMPASSContainer(data, counts, meta, stimulation_id,
    individual_id, sample_id)
}
\arguments{
  \item{data}{A list of matrices. Each matrix \code{M_i} is
  made up of \code{N_i} cells by \code{K} markers; for
  example, it could be the intensity information from an
  intracellular cytokine experiment.  Each element of the
  list should be named; this name denotes which sample the
  cell intensities were measured from.}

  \item{counts}{A named integer vector of the cell counts
  for each sample in \code{data}. This element is required
  if 'null' cells; i.e., cells that expressed no markers,
  are not included in \code{data}.}

  \item{meta}{A \code{data.frame} of metadata, describing
  the individuals in the experiment. Each row in
  \code{meta} should correspond to a row in \code{data}.
  There should be one row for each sample; i.e., one row
  for each element of \code{data}.}

  \item{stimulation_id}{The name of the vector in
  \code{meta} that denotes the type of stimulation each
  sample received.}

  \item{individual_id}{The name of the vector in
  \code{meta} that denotes the individuals from which
  samples were drawn.}

  \item{sample_id}{The name of the vector in \code{meta}
  that denotes the samples.  This vector should contain all
  of the names in the \code{data} input.}
}
\description{
  This function generates the data container suitable for
  use with \code{COMPASS}.
}

