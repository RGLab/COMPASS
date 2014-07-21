#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector boost_digamma(NumericVector x) {
  int n = x.size();
  NumericVector output = no_init(n);
  for (int i=0; i < n; ++i) {
    output[i] = boost::math::digamma(x[i]);
  }
  return output;
}

// [[Rcpp::export]]
NumericVector R_digamma(NumericVector x) {
  return digamma(x);
}

/* DIGAMMA.C - Compute the digamma function. */

/* Copyright (c) 1995-2004 by Radford M. Neal
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


/* digamma(x) is defined as (d/dx) log Gamma(x).  It is computed here
   using an asymptotic expansion when x>5.  For x<=5, the recurrence
   relation digamma(x) = digamma(x+1) - 1/x is used repeatedly.  See
   Venables & Ripley, Modern Applied Statistics with S-Plus, pp. 151-152. */


/* COMPUTE THE DIGAMMA FUNCTION.  Returns -inf if the argument is an integer
   less than or equal to zero. */

double digamma_radford_neal (double x)
{
  double r, f, t;

  r = 0;

  while (x<=5)
  { r -= 1/x;
    x += 1;
  }

  f = 1/(x*x);

  // t = f*(-1/12.0 + f*(1/120.0 + f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0)))));
  t = f*(-1/12.0 + f*(1/120.0 + f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0
       + f*(691/32760.0 + f*(-1/12.0 + f*3617/8160.0)))))));

  return r + log(x) - 0.5/x + t;
}

double digamma_mod(double x) {
  double r, f1, f2, f3, f4;
  r = 0;

  while (x < 5) {
    r -= 1/x;
    x += 1;
  }

  f1 = 1 / (x*x);
  f2 = f1 * f1;
  f3 = f2 * f1;
  f4 = f2 * f2;

  static const double c1 = -1/12.0;
  static const double c2 = 1/120.0;
  static const double c3 = -1/252.0;
  static const double c4 = 1/240.0;
  static const double c5 = -1/132.0;
  static const double c6 = 691/32760.0;
  static const double c7 = -1/12.0;
  static const double c8 = 3617/8160.0;

  double A = c1 + c2*f1 + c3*f2 + c4*f3;
  double B = c5*f1 + c6*f2 + c7*f3 + c8*f4;
  return r + log(x) - 0.5 / x + f1*A + f4*B;
}

// [[Rcpp::export]]
NumericVector c_digamma(NumericVector x) {
  int n = x.size();
  NumericVector output = no_init(n);
  for (int i=0; i < n; ++i) {
    output[i] = digamma_radford_neal(x[i]);
  }
  return output;
}

// [[Rcpp::export]]
NumericVector c_digamma_mod(NumericVector x) {
  int n = x.size();
  NumericVector output = no_init(n);
  for (int i=0; i < n; ++i) {
    output[i] = digamma_mod(x[i]);
  }
  return output;
}

/*** R
library(microbenchmark)
x <- rgamma(1E4, 1E-1)
microbenchmark(
  R_digamma(x),
  boost_digamma(x),
  c_digamma(x),
  c_digamma_mod(x)
)
x <- rgamma(1E4, 10)
microbenchmark(
  R_digamma(x),
  boost_digamma(x),
  c_digamma(x),
  c_digamma_mod(x)
)
x <- seq(1E-5, 1, by=1E-4)
R <- R_digamma(x)
B <- boost_digamma(x)
C <- c_digamma(x)
D <- c_digamma_mod(x)
mean(R - B)
mean(R - C)
mean(R - D)
max( R - C )
all.equal(R, C)
x <- rgamma(1E4, 1E-1)
R <- R_digamma(x)
B <- boost_digamma(x)
C <- c_digamma(x)
D <- c_digamma_mod(x)
mean(R - B)
mean(R - C)
mean(R - D)
max( R - C )
all.equal(R, C)
all.equal(R, D)
*/
