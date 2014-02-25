#include <Rcpp.h>
using std::vector;

void updatealpha_scalar(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM, int K1, int t1, 
                        double xp_var, double sqrt_var1, double sqrt_var2, int xtt, vector<int>& xgammat, Rcpp::IntegerMatrix& xd, 
                        Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, Rcpp::NumericMatrix& pt1, Rcpp::NumericMatrix& pt2, Rcpp::NumericMatrix& pt3,
                        double xlambda, double xbeta, Rcpp::NumericVector& alpha, double xsig_alpha1, double xalpha1, Rcpp::IntegerVector& xAalpha);
