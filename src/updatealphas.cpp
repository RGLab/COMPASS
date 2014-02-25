#include <Rcpp.h>

extern "C" double digamma(double);

using namespace std;
using namespace Rcpp;


void updatealphas(vector<double>& xalphast,vector<int>& xn_s, int xK, int xI, vector<double>& xlambda_s, vector<int>& xgammat, 
                  Rcpp::NumericVector& sqrt_var1,Rcpp::NumericVector& sqrt_var2, vector<double>& xp_var, int xtt, vector<int>& xAalphas)
{

    double delF = 0.0;
    double psik = 0.0;
    double log1 = 0.0;
    double log2 = 0.0;  int flag1=0;   int flagkk=0; int lp1=0;
    for (int kk = 0; kk < xK; kk++) {
        delF = 0.0;  
        psik = digamma(xalphast[kk]);
        log1 = 0.0;
        log2 = 0.0;
        for (int i = 0; i < xI; i++) {
            lp1=0;
            for (int k = 0; k < xK; k++) {
              if (xgammat[i+xI*k] == 1) { lp1 +=1;}
            }
            vector<int> p1(lp1); flag1=0; flagkk=0;
            for (int k = 0; k < xK; k++) {
              if (xgammat[i+xI*k] == 1) { 
                p1[flag1] = k;
                flag1 += 1;
                if (k == kk) {flagkk = 1;} 
              }       
            }
            double sum_alp_ns = 0.0;
            double sum_alp = 0.0;
            double sum_gl_alp = 0.0;
            double sum_gl_alp_ns = 0.0;
            for (int k = 0; k<lp1; k++){
                double sums = xalphast[p1[k]] + xn_s[i+xI*p1[k]];
                sum_alp_ns += sums;
                sum_alp += xalphast[p1[k]];
                sum_gl_alp += lgamma(xalphast[p1[k]]);
                sum_gl_alp_ns += lgamma(sums);
            }
            if (flagkk > 0) {
               delF += digamma(xn_s[i+xI*kk]+xalphast[kk]) - psik - digamma(sum_alp_ns) + digamma(sum_alp);
            }
            if (lp1>0) {
               log2 += -(sum_gl_alp-lgamma(sum_alp)) + (sum_gl_alp_ns - lgamma(sum_alp_ns));
            }
        }
        double mean_p = std::max(0.01, xalphast[kk]+delF/xtt);
        Rcpp::NumericVector alpha_s_p= Rcpp::rnorm(1, mean_p, sqrt_var1[kk]);
        if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var[kk])) == 1) {
            alpha_s_p= Rcpp::rnorm(1, mean_p, sqrt_var1[kk]);}
        else { alpha_s_p = Rcpp::rnorm(1, mean_p, sqrt_var2[kk]);}
        
        if (alpha_s_p[0]>0.0 && alpha_s_p[0]<=xlambda_s[kk]) {
           vector<double> alp(xK);
        
           for (int i = 0; i<xK; i++) {
               alp[i] = xalphast[i];
           }
           alp[kk] = alpha_s_p[0];
           // log2 += log(xp_var[kk]*gsl_ran_gaussian_pdf(alp[kk]-mean_p, sqrt_var1[kk])+(1-xp_var[kk])*gsl_ran_gaussian_pdf(alp[kk]-mean_p, sqrt_var2[kk])); 
           log2 += log(
             xp_var[kk] * Rf_dnorm4(alp[kk], mean_p, sqrt_var1[kk], 0) +
             (1 - xp_var[kk]) * Rf_dnorm4(alp[kk], mean_p, sqrt_var2[kk], 0)
           );
           
           delF = 0.0; psik = digamma(alp[kk]);
           for ( int i = 0; i < xI; i++) {
               
              lp1 = 0; 
              for (int k = 0; k < xK; k++) {
                  if (xgammat[i+xI*k] == 1) { lp1 +=1;}       
               }

              vector<int> p1(lp1); flag1 =0; flagkk = 0;
               for (int k = 0; k < xK; k++) {
                  if (xgammat[i+xI*k] == 1) { 
                    p1[flag1] = k;
                     flag1 += 1;
                     if (k == kk) {flagkk = 1;}
                  }       
               }

               double sum_alp_ns = 0.0;
               double sum_alp = 0.0;
               double sum_gl_alp = 0.0;
               double sum_gl_alp_ns = 0.0;
               for (int k = 0; k<lp1; k++){
                   double sums = alp[p1[k]] + xn_s[i+xI*p1[k]];
                   sum_alp_ns += sums;
                   sum_alp += alp[p1[k]];
                   sum_gl_alp += lgamma(alp[p1[k]]);
                   sum_gl_alp_ns += lgamma(sums);
               }
               if (flagkk > 0) {
                   delF += digamma(xn_s[i+xI*kk]+alp[kk]) - psik - digamma(sum_alp_ns) + digamma(sum_alp);
               }
               if (lp1>0) {
                   log1 += -(sum_gl_alp-lgamma(sum_alp)) + (sum_gl_alp_ns - lgamma(sum_alp_ns));
               }
           }
           mean_p = std::max(0.01, alp[kk] + delF/xtt);
           // log1 +=log(xp_var[kk]*gsl_ran_gaussian_pdf(xalphast[kk]-mean_p, sqrt_var1[kk])+(1-xp_var[kk])*gsl_ran_gaussian_pdf(xalphast[kk]-mean_p, sqrt_var2[kk])); 
           log1 += log(
             xp_var[kk] * Rf_dnorm4(xalphast[kk], mean_p, sqrt_var1[kk], 0) +
             (1 - xp_var[kk]) * Rf_dnorm4(xalphast[kk], mean_p, sqrt_var2[kk], 0)
           );
             
           if (log(Rcpp::as<double>(Rcpp::runif(1))) <= (log1 - log2)) {
                xalphast[kk] = alp[kk];
                xAalphas[kk] = 1;
           } 
        }
    
    }

}
