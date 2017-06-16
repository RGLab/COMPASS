// [[Rcpp::interfaces(r)]]

#include <Rcpp.h>
using namespace Rcpp;

extern "C" SEXP mat2vec(SEXP);

// [[Rcpp::export]]
IntegerMatrix CellCounts_character(List data, List combinations) {

  Function list2env("list2env");

  int m = data.size();          // rows
  int n = combinations.size();  // columns
  IntegerMatrix output(m, n);

  // loop over data
  for (int i = 0; i < m; ++i) {

    SEXP m = PROTECT(mat2vec(as<LogicalMatrix>(data[i])));
    Environment env = list2env(m);

    // loop over combinations
    for (int j = 0; j < n; ++j) {
      ExpressionVector expr = as<ExpressionVector>(combinations[j]);
      SEXP result = PROTECT(Rf_eval(expr[0], env));
      int* result_ptr = INTEGER(result);
      for (int k = 0; k < Rf_length(result); ++k) {
        output(i, j) += result_ptr[k];
      }
      UNPROTECT(1);
    }

    UNPROTECT(1);
  }

  return output;
}
