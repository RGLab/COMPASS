##' Fit the discrete COMPASS Model
##' 
##' This function fits the \code{COMPASS} model from a user-provided set of
##' stimulated / unstimulated matrices.
##' 
##' @param n_s The cell counts for stimulated cells.
##' @param n_u The cell counts for unstimulated cells.
##' @param meta A \code{data.frame} of metadata, describing the individuals
##'   in the experiment. Each row in \code{meta} should correspond to a row
##'   in \code{data}. There should be one row for each sample;
##'   i.e., one row for each element of \code{n_s} and \code{n_u}.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn. 
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @param iterations The number of iterations (per 'replication') to perform.
##' @param replications The number of 'replications' to perform. In order to
##'   conserve memory, we only keep the model estimates from the last replication.
##' @param verbose Boolean; if \code{TRUE} we output progress information.
##' @return A \code{list} with class \code{COMPASSResult} with two components,
##'   the \code{fit} containing parameter estimates and parameter acceptance
##'   rates, and \code{data} containing the generated data used as input for
##'   the model.
##' @examples \dontrun{
##' set.seed(123)
##' n <- 10 ## number of samples
##' k <- 3 ## number of markers
##' 
##' ## generate some sample data
##' sid_vec <- paste0("sid_", 1:n) ## sample ids; unique names used to denote samples
##' iid_vec <- rep_len( paste0("iid_", 1:(n/2) ), n ) ## individual ids
##' data <- replicate(n, {
##'   nrow <- round(runif(1) * 1E4 + 1000)
##'   ncol <- k
##'   vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
##'   vals[ vals < 2000 ] <- 0
##'   output <- matrix(vals, nrow, ncol)
##'   output <- output[ apply(output, 1, sum) > 0, ]
##'   colnames(output) <- paste0("M", 1:k)
##'   return(output)
##' })
##' meta <- data.frame(
##'   sid=sid_vec,
##'   iid=iid_vec,
##'   trt=rep( c("Control", "Treatment"), each=(n/2) )
##' )
##' 
##' ## generate counts for n_s, n_u
##' n_s <- CellCounts( data[1:(n/2)], Combinations(k) )
##' n_u <- CellCounts( data[(n/2+1):n], Combinations(k) )
##' 
##' ## A smaller number of iterations is used here for running speed;
##' ## prefer using more iterations for a real fit
##' SimpleCOMPASS(n_s, n_u, meta, "iid", "sid", iterations=100)
##' }
SimpleCOMPASS <- function(n_s, n_u, meta, individual_id, sample_id,
  iterations=1E4, replications=8, verbose=TRUE) {
  
  if (!all(colnames(n_s) == colnames(n_u))) {
    stop("The column names of 'n_s' and 'n_u' do not match.")
  }
  
  n_markers <- log2( ncol(n_s) )
  if (!(n_markers == as.integer(n_markers))) {
    stop("Could not infer the number of markers correctly; be sure that ",
      "each possible combination of cells is represented in your counts matrix.")
  }
  
  ## Guess the marker names
  marker_names <- unique(
    unlist( strsplit( gsub("!", "", colnames(n_s)), "&", fixed=TRUE ) )
  )
  
  if (length(marker_names) != n_markers) {
    stop("Internal error: could not infer the marker names from the ",
      "data supplied!", call.=FALSE)
  }
  
  cats <- as.data.frame( matrix(0, nrow=ncol(n_s), ncol=n_markers) )
  rownames(cats) <- colnames(n_s)
  colnames(cats) <- 
  for (i in seq_along(cats)) {
    cats[, i] <- as.integer(grepl( paste0( colnames(cats)[i], "+" ), rownames(cats), fixed=TRUE ))
  }
  cats$Counts <- apply(cats, 1, sum)
  cats <- as.matrix(cats)
  
  n_s <- as.matrix(n_s)
  counts_s <- rowSums(n_s)
  
  n_u <- as.matrix(n_u)
  counts_u <- rowSums(n_u)
  
  .fit <- .COMPASS.discrete(
    n_s=n_s,
    n_u=n_u,
    categories=cats,
    iterations=100,
    replications=8,
    verbose=TRUE
  )
  
  fit <- list(
    fit=.fit,
    data=list(
      n_s=n_s,
      n_u=n_u,
      counts_s=counts_s,
      counts_u=counts_u,
      categories=cats,
      meta=meta,
      sample_id=sample_id,
      individual_id=individual_id
    )
  )
  
  class(fit) <- c("SimpleCOMPASSResult", "COMPASSResult")
  return(fit)
  
}
