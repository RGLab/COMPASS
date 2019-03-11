#' Compute a response probability from COMPASS mcmc samples.
#'
#' @param x  a \code{COMPASSResult} object.
#' @param markers a \code{vector} of marker names.
#' @param degree the \code{numeric} degree of functionality to test.
#' @description
#' Compute a response probability based on the selected markers, evaluating the probability
#' that a subject exhibits a response of size \code{degree} or greater.
#' i.e., the probability of at least \code{degree} markers exhibiting an antigen specific response.
#'
#' @details
#' The response is computed from the sampled Gamma matrix, subsetting on the selected markers, and
#' @return  A \code{vector} of response probabilities for each subject.
#' @export
#'
#' @examples
#' Response(CR, markers = c("M1","M2","M3"), degree = 2)
Response <- function(x, markers, degree){
  UseMethod("Response")
}

##' @rdname FunctionalityScore
##' @export
Response.COMPASSResult <- function(x, markers = NULL, degree = 1) {
  ## we drop the last column as it is the 'NULL' category
  if (is.null(markers)) {
    markers <- markers(x)
  }
  if (degree > length(markers)) {
    stop("Invalid degree: ", degree, ". Only ", length(markers)," markers provided")
  }
    if (!all(markers %in% markers(x))) {
      stop("Invalid marker names")
    }
    message("Computing the probability of response of degree >= ",degree, " from markers: ", paste(markers, collapse = ", "))
    new_categories = unique(categories(x, FALSE)[, markers, drop = FALSE])
    all_categories = categories(x, FALSE)[, markers, drop = FALSE]
    suppressWarnings({dmat = as.matrix(pdist(new_categories, all_categories))})
    cat_indices = apply(dmat, 1, function(y)
      which(y == 0))
    if (!is.matrix(cat_indices)) {
      cat_indices <- matrix(cat_indices, ncol = length(cat_indices))
    }
    new_mean_gamma = apply(cat_indices, 2, function(i)
      apply(Gamma(x)[, i, ], 1, mean))
    new_categories = cbind(new_categories, Counts = rowSums(new_categories))
    reord = c(setdiff(1:nrow(new_categories), which(new_categories[, "Counts"] ==
                                                      0)), which(new_categories[, "Counts"] == 0))
    new_categories = new_categories[reord, ]
    new_mean_gamma = new_mean_gamma[, reord]
    colnames(new_mean_gamma) = apply(new_categories[, -ncol(new_categories)], 1, function(x)
      paste0(x, collapse = ""))
    include_cols <- new_categories[,"Counts"] >= degree
    if (sum(include_cols) == 0) {
      response <- matrix(rep(0,nrow(new_mean_gamma)),nrow=nrow(new_mean_gamma),ncol=1)
      rownames(response) <- rownames(new_mean_gamma)
      colnames(response) <- paste0("Pr(response|degree >=",degree,")")
    }else{
      response <- new_mean_gamma[,include_cols, drop = FALSE]
      response <- rowMeans(response)
      response <- matrix(response, ncol = 1, nrow = length(response))
      rownames(response) <- rownames(new_mean_gamma)
      colnames(response) <- paste0("Pr(response | degree >=",degree,")")
    }
    response
}

