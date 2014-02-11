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

compute_posterior <- function(x, as.matrix=FALSE) {
  
  output <- lapply( 1:nrow(x$data$n_s), function(i) {
    .Call( C_samplePuPs,
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

##' @rdname Posterior
##' @export
PosteriorDiff <- function(x) {
  
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' must be an object of class 'COMPASSResult'")
  }
  
  post <- x$fit$posterior
  output <- sapply(post, "[[", "diff")
  nm <- colnames( x$data$n_s )
  rownames(output) <- nm[ -length(nm) ]
  return( t(output) )
  
}

##' @rdname Posterior
##' @export
PosteriorLogDiff <- function(x) {
  
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' must be an object of class 'COMPASSResult'")
  }
  
  post <- x$fit$posterior
  output <- sapply(post, "[[", "logd")
  nm <- colnames( x$data$n_s )
  rownames(output) <- nm[ -length(nm) ]
  return( t(output) )
  
}

##' @rdname Posterior
##' @export
PosteriorPs <- function(x) {
  
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' must be an object of class 'COMPASSResult'")
  }
  
  post <- x$fit$posterior
  output <- sapply(post, "[[", "p_s")
  nm <- colnames( x$data$n_s )
  rownames(output) <- nm[ -length(nm) ]
  return( t(output) )
  
}

##' @rdname Posterior
##' @export
PosteriorPu <- function(x) {
  
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' must be an object of class 'COMPASSResult'")
  }
  
  post <- x$fit$posterior
  output <- sapply(post, "[[", "p_u")
  nm <- colnames( x$data$n_s )
  rownames(output) <- nm[ -length(nm) ]
  return( t(output) )
  
}
