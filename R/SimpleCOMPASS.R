##' Fit the discrete COMPASS Model
##'
##' This function fits the \code{COMPASS} model from a user-provided set of
##' stimulated / unstimulated matrices. See the NOTE for important details.
##'
##' @note n_s and n_u counts matrices should contain ALL 2^M possible combinations of markers, even if they are 0 for some combinations. The code expects the marker combinations to be named in the following way:
##' \code{"M1&M2&!M3"} means the combination represents cells expressing marker "M1" and "M2" and not "M3". For 3 markers, there should be 8 such combinations, such that n_s and n_u have 8 columns.
##' @param n_s The cell counts for stimulated cells.
##' @param n_u The cell counts for unstimulated cells.
##' @param meta A \code{data.frame} of metadata, describing the individuals
##'   in the experiment. Each row in \code{meta} should correspond to a row
##'   in \code{data}. There should be one row for each subject;
##'   i.e., one row for each element of \code{n_s} and \code{n_u}.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn.
##' @param iterations The number of iterations (per 'replication') to perform.
##' @param replications The number of 'replications' to perform. In order to
##'   conserve memory, we only keep the model estimates from the last replication.
##' @param verbose Boolean; if \code{TRUE} we output progress information.
##' @param seed A seed for the random number generator. Defaults to 100.
##' @return A \code{list} with class \code{COMPASSResult} with two components,
##'   the \code{fit} containing parameter estimates and parameter acceptance
##'   rates, and \code{data} containing the generated data used as input for
##'   the model.
##' @export
##' @examples
##'  set.seed(123)
##' n <- 10 ## number of subjects
##' k <- 3 ## number of markers
##'
##' ## generate some sample data
##' iid_vec <- paste0("iid_", 1:n) # Subject id
##' data <- replicate(2*n, {
##' nrow <- round(runif(1) * 1E4 + 1000)
##' ncol <- k
##' vals <- rexp( nrow * ncol, runif(1, 1E-5, 1E-3) )
##' vals[ vals < 2000 ] <- 0
##' output <- matrix(vals, nrow, ncol)
##'output <- output[ apply(output, 1, sum) > 0, ]
##'colnames(output) <- paste0("M", 1:k)
##'return(output)
##'})
##'
##' meta <- cbind(iid=iid_vec, data.frame(trt=rep( c("Control", "Treatment"), each=n/2 )))
##'
##' ## generate counts for n_s, n_u
##' n_s <- CellCounts( data[1:n], Combinations(k) )
##' n_u <- CellCounts( data[(n+1):(2*n)], Combinations(k) )
##' rownames(n_s) = unique(meta$iid)
##' rownames(n_u) = rownames(n_s)

##' ## A smaller number of iterations is used here for running speed;
##' ## prefer using more iterations for a real fit
##' scr = SimpleCOMPASS(n_s, n_u, meta, "iid", iterations=1000)
SimpleCOMPASS <- function(n_s, n_u, meta, individual_id,
  iterations=1E4, replications=8, verbose=TRUE,seed=100) {


  # Order, n_s, n_u, and meta (if needed)
  rn_s <- rownames(n_s)
  rn_u <- rownames(n_u)
  iid <- as.character(meta[, individual_id])

  if(!identical(rn_s, rn_u) | !identical(rn_s, iid)) {
  n_s <- n_s[order(rn_s),]
  n_u <- n_u[order(rn_u),]
  meta <- meta[order(iid),]
  message("Ordering meta, n_s and n_u by individual_id since this wasn't done.\n",
          "If you think this is an error, check your data and rerun the code.")
  }
  cat("Setting the seed to ",seed,"\n");
  set.seed(seed)
  if (!all(colnames(n_s) == colnames(n_u))) {
    stop("The column names of 'n_s' and 'n_u' do not match.")
  }

  n_markers <- log2( ncol(n_s) )
  if (!(n_markers == as.integer(n_markers))) {
    warning("Could not infer the number of markers correctly; it looks like ",
      "you may have filtered some cell-subsets. If that is the case, you can ignore this warning.")
  }

  ## Guess the marker names
  marker_names <- unique(
    unlist( strsplit( gsub("!", "", colnames(n_s)), "&", fixed=TRUE ) )
  )
  n_markers <- length(marker_names)
  cats <- as.data.frame( matrix(0, nrow=ncol(n_s), ncol=n_markers) )
  rownames(cats) <- colnames(n_s)
  colnames(cats) = marker_names

  for (i in seq_along(cats)) {
    #cats[, i] <- as.integer(grepl( paste0( colnames(cats)[i], "+" ), rownames(cats), fixed=TRUE ))
    cats[,i] <-
      as.integer(!grepl(paste0("!",colnames(cats)[i],"(&|$)+"),rownames(cats),fixed =
                          FALSE))
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
    iterations=iterations,
    replications=replications,
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
      individual_id=individual_id
    )
  )
  class(fit) <- c("COMPASSResult")
  return(fit)
}
