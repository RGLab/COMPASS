#' Compute a response probability from COMPASS mcmc samples.
#'
#' @param x  a \code{COMPASSResult} object.
#' @param markers a \code{vector} of marker names.
#' @param degree the \code{numeric} degree of functionality to test.
#' @param max.prob \code{logical} Use the max probability rather than the average across subsets. Defaults to FALSE.
#' @param cond_response \code{logical} Renormalize the probabilities to condition on response. i.e. given that there
#' is a response, what's the probability of a response of degree X or greater?
#' @description
#' Compute a response probability based on the selected markers, evaluating the probability
#' that a subject exhibits a response of size \code{degree} or greater.
#' i.e., the probability of at least \code{degree} markers exhibiting an antigen specific response.
#'
#' @details
#' The response is computed from the sampled Gamma matrix, subsetting on the selected markers, and
#' @return  A \code{vector} of response probabilities for each subject.
#' @examples
#'
#' Response(CR, markers = c("M1","M2","M3"), degree = 2)
#' @export
Response <- function(x, markers, degree, max.prob, cond_response){
  UseMethod("Response")
}

##' @rdname Response
##' @export
Response.COMPASSResult <- function(x, markers = NULL, degree = 1, max.prob = FALSE, cond_response = FALSE) {
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
    # new_gamma <- Gamma(x)[,cat_indices,]
    new_categories = cbind(new_categories, Counts = rowSums(new_categories))
    reord = c(setdiff(1:nrow(new_categories), which(new_categories[, "Counts"] ==
                                                      0)), which(new_categories[, "Counts"] == 0))
    new_categories = new_categories[reord, ]
    new_mean_gamma = new_mean_gamma[, reord]
    # new_gamma = new_gamma[,reord,]
    colnames(new_mean_gamma) = apply(new_categories[, -ncol(new_categories)], 1, function(x)
      paste0(x, collapse = ""))
    # dimnames(new_gamma)[[2]] = apply(new_categories[, -ncol(new_categories)], 1, function(x)
      # paste0(x, collapse = ""))
    include_cols <- new_categories[,"Counts"] >= degree
    if (sum(include_cols) == 0) {
      response <- matrix(rep(0,nrow(new_mean_gamma)),nrow = nrow(new_mean_gamma),ncol=1)
      rownames(response) <- rownames(new_mean_gamma)
      colnames(response) <- paste0("Pr(response|degree >=",degree,")")
    }else{
      #condition on degree >= x
      response <- new_mean_gamma[,include_cols, drop = FALSE]
      if (!max.prob & !cond_response) {
        # probability of response from total probability
        # response <- structure(
        #   rowMeans(new_gamma[, include_cols, , drop = FALSE]) / (
        #   rowMeans(new_gamma[, !include_cols, , drop = FALSE]) +
        #     rowMeans(new_gamma[, include_cols, , drop = FALSE]) +
        #     rowMeans(1 - new_gamma[, !include_cols, , drop = FALSE]) +
        #     rowMeans(1 - new_gamma[,include_cols, , drop = FALSE])),
        #   dim = c(nrow(response),1),
        #   names = NULL,
        #   dimnames = list(rownames(response),paste0("Pr(response | degree >=",degree,")")))

        response <- structure(rowMeans(response), dim = c(nrow(response),1), names = NULL, dimnames = list(rownames(response),paste0("Pr(response | degree >=",degree,")")))
      }else if (max.prob & !cond_response){
        response <- structure(apply(response,1,max),dim=c(nrow(response),1), names = NULL, dimnames = list(rownames(response), paste0("Pr(response | degree >=",degree,")")))
      }else if (cond_response & ! max.prob){
        response <-
          structure(
            rowSums(exp(log(
              apply(Gamma(x), 1:2, mean)
            ) - log(
              rowSums(apply(Gamma(x), 1:2, mean)) + 1e-10
            ))[, include_cols, drop = FALSE]),
            dim = c(nrow(response), 1),
            names = NULL,
            dimnames = list(
              rownames(response),
              paste0("Pr(degree >= ",degree," | response )")
            )
          )
      }
      # response <- matrix(response, ncol = 1, nrow = length(response))
      # rownames(response) <- rownames(new_mean_gamma)
      # colnames(response) <- paste0("Pr(response | degree >=",degree,")")
    }
    response
}

