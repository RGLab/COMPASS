// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;

// [[Rcpp::export]]
double Rcpp_add_gaussian(NumericVector x, double mean, double sd) {
  NumericVector dnorms = dnorm(x, mean, sd, true);
  double output = sum(dnorms);
  return output;
}

// [[Rcpp::export]]
double gsl_add_gaussian(NumericVector x, double mean, double sd) {
  double sum = 0;
  for (int i=0; i < x.size(); ++i) {
    sum += gsl_ran_gaussian_pdf(x[i] - mean, sd);
  }
  return sum;
}

// [[Rcpp::export]]
List gsl_compare(double x, double mean, double sd) {
  double num1 = gsl_ran_gaussian_pdf(x - mean, sd);
  NumericVector x_(1);
  x_[0] = x;
  NumericVector num2 = dnorm(x_, mean, sd);
  
  double num3 = Rf_dnorm4(x, mean, sd, 0);
  double num4 = Rf_dnorm4(x, mean, sd, 1);
  double num5 = log(num3);
  return List::create( _["GSL"] = num1, _["Rcpp"] = num2, _["R"] = num3,
    _["R_log"] = num4, _["log[R]"] = num5);
}

// [[Rcpp::export]]
NumericVector gsl_exp(NumericVector x, double rate) {
  NumericVector output(x.size());
  for (int i=0; i < x.size(); ++i) {
    output[i] = gsl_ran_exponential_pdf( x[i], rate );
  }
  return output;
}

// [[Rcpp::export]]
NumericVector R_exp(NumericVector x, double rate) {
  return dexp(x, rate);
}

// [[Rcpp::export]]
List gsl_compare_exp(double x, double rate) {
  
  return List::create(
    _["GSL"] = gsl_ran_exponential_pdf(x, rate),
    _["R"] = Rf_dexp(x, rate, 0)
  );
  
}

/*** R
library(microbenchmark)
x <- 1
mu <- 2
sd <- 1
gsl_compare(x, mu, sd)

x <- rnorm(1E5)
Rcpp_add_gaussian(x, mu, sd)
microbenchmark(
  Rcpp_add_gaussian(x, mu, sd),
  gsl_add_gaussian(x, mu, sd)
)
*/
