##' Compute the Posterior Difference, Log Ratio
##' 
##' Computes the posterior difference and log ratio from a \code{COMPASSResult}.
##' 
##' @param x An object of class \code{COMPASSResult}.
##' @export
Posterior <- function(x) {  
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' must be an object of class 'COMPASSResult'")
  }
  
  return(x$fit$posterior)
}

compute_posterior <- function(x) {
  output <- lapply( 1:nrow(x$data$n_s), function(i) {
    .Call( "samplePuPs",
      x$fit$alpha_u,
      x$fit$alpha_s,
      x$fit$gamma[i, , ],
      dim(x$fit$gamma)[[3]],
      dim(x$fit$gamma)[[2]],
      x$data$n_s[i, ],
      x$data$n_u[i, ],
      x$fit$categories,
      ncol(x$fit$categories) - 1L
    )
  })
  
  names(output) <- rownames(x$data$n_s)
  return(output)
}
