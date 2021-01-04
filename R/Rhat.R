#'@title Diagnostic of a set of COMPASS Models.
#' @param x a list of compass model fits of the same data with the same number of iterations, different seeds.
#' Run some mcmc diagnostics on a series of COMPASS model fits.
#' Assuming the input is a list of model fits for the same data with the same number of iterations and different seeds.
#' Run Gelman's Rhat diagnostics on the alpha_s and alpha_u hyperparameter chains, treating each model as an independent chain.
#' Rhat should be near 1 but rarely are in practice. Very large values may be a concern.
#' The method returns an average model, by averaging the mean_gamma matrices (equally weighted since each input has the same number of iterations).
#' This mean model should be better then any of the individual models.
#' It can be plotted via "plot(result$mean_model)".
#' @importFrom stats fisher.test
#' @export
COMPASSMCMCDiagnosis<-function(x){
    diag<-list()
    diag$alpha_s<-coda::gelman.diag(Map(function(x)coda::as.mcmc(x$fit$alpha_s),x))
    diag$alpha_u<-coda::gelman.diag(Map(function(x)coda::as.mcmc(x$fit$alpha_u),x))
  mean_gamma <- apply(Map(function(x) abind(x, along = 3), Map(function(x) x$fit$mean_gamma, x))[[1]], 1:2, mean)
  mean_model <- x[[1]]
  mean_model$fit$mean_gamma <- mean_gamma
    return(list(diag=diag,mean_model=mean_model))
}

