#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <time.h>

#include <R_ext/Utils.h>


using std::vector;
void updatealphau_RW(vector<double>& xalphaut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, vector<double>& xlambda_u, Rcpp::NumericVector& sqrt_var, int xtt, vector<int>& xgammat,vector<int>& xAalphau);
