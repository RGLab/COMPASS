#include <Rcpp.h>
#include <stdio.h>

// cpp codes for generate posterior ps and pu for each subject i
RcppExport SEXP samplePuPs(SEXP alphau, SEXP alphas, SEXP gammat,  SEXP T, SEXP K , SEXP nsi, SEXP nui, SEXP d, SEXP M)
{
    BEGIN_RCPP
    
    Rcpp::IntegerVector xnsi(nsi); // n_s[i,]
    Rcpp::IntegerVector xnui(nui); //n_u[i,]
    Rcpp::NumericMatrix xalphau(alphau); // posterior samples for alpha_u
    Rcpp::NumericMatrix xalphas(alphas); // posterior samples for alpha_s
    Rcpp::IntegerMatrix xgammat(gammat); // posterior samples for gamma_i
    Rcpp::IntegerMatrix xd(d); //cytokine combination indicator matrix
    
    int xM = Rcpp::as<int>(M); // # markers
    int xT = Rcpp::as<int>(T); // # MCMC iterations used
    int xK = Rcpp::as<int>(K); // # Categories
    int K1 = xK-1; 
    Rcpp::NumericVector xpu(xK);
    Rcpp::NumericVector xps(xK);

    Rcpp::NumericVector DIFF(K1); // Ps-Pu
    Rcpp::NumericVector LogD(K1); // LogPs-logPu
      
  
    int sum = 0;
    double alpha_u[xK];
    double alpha_s[xK];
    
    Rcpp::NumericMatrix psi(xT,xK);
    Rcpp::NumericMatrix pui(xT,xK);
    double diff[xK];
    double logd[xK];
    double tmp = 0.;
    Rcpp::RNGScope scope;
    double sau = 0.;
    double sas = 0.;
    double sass = 0.;

    for (int tt = 0; tt < xT; tt++) {
      sau = 0.;
      sas = 0.;
      for ( int j=0; j<xK; j++){
        alpha_u[j] = xalphau(tt,j) + xnui[j];
        alpha_s[j] = xalphas(tt,j) + xnsi[j];
        sau += alpha_u[j];
        sas += alpha_s[j];
      }   
      sum = 0;
      for (int j = 0; j < xK; j++) {
          pui(tt,j) = alpha_u[j]/sau;
          sum += xgammat(j,tt);
      }
      if (sum==0){
        for(int j=0; j<xK; j++){
          psi(tt,j) = pui(tt,j);
          diff[j] = 0;
          logd[j] = 0;
         
        }
      } else if (sum==xK) {
        for(int j=0; j<xK; j++) {
          psi(tt,j) = alpha_s[j]/sas;
          diff[j] = psi(tt,j) - pui(tt,j);
          logd[j] = log2(psi(tt,j)) - log2(pui(tt,j));
         
        }
      }else{
        int l0 = xK-sum;
        int place0[l0];   
        int place1[sum];
        int flag0 = 0; int flag1 = 0;
        for ( int j=0; j<xK; j++){
          if (xgammat(j,tt) == 0) {
              place0[flag0] = j;
              flag0 += 1;
          } else {
              place1[flag1] = j;
              flag1 += 1;
          }   
        }
        tmp = 0.;
        for(int j = 0; j<l0; j++){
          psi(tt,place0[j]) = pui(tt,place0[j]);
          tmp += psi(tt,place0[j]);
          diff[place0[j]] = 0;
          logd[place0[j]] = 0;
          
        }
        double as[sum];
        sass = 0.;
        for (int j=0; j<sum;j++){
       
          as[j] = alpha_s[place1[j]];
          sass += as[j];
        }

        for(int j=0; j<sum; j++){
          psi(tt,place1[j]) = as[j]*(1-tmp)/sass;
          diff[place1[j]] = psi(tt,place1[j]) - pui(tt,place1[j]);
          logd[place1[j]] = log2(psi(tt,place1[j])) - log2(pui(tt,place1[j]));
        }
       
      }
       for(int j=0; j<K1; j++){
         LogD[j] += logd[j]/xT;
         DIFF[j] += diff[j]/xT;   
       }

    }
   
    return Rcpp::List::create(Rcpp::Named("diff") = DIFF, Rcpp::Named("logd") = LogD);
    END_RCPP


}



