#include <Rcpp.h>
using namespace Rcpp;

#define rownames(x) as<CharacterVector>(VECTOR_ELT( x.attr("dimnames"), 0 ))
#define colnames(x) as<CharacterVector>(VECTOR_ELT( x.attr("dimnames"), 1 ))

// x is a list of logical matrices, corresponding to 'expressed' or 'unexpressed'
// with rows as cells and columns as cytokines

// combos is a list of integer vectors as follows:
// for each cytokine we wish to check, we pass its index (1, 2, ...)
// if we want to check if it is positive, we pass a positive number
// if we want to check if it is negative, we pass a negative number, of the
// same index

// eg: for a matrix with columns CD4 and TNFa, we might pass combos=list(c(1, -2))
// to get counts that are CD4+TNFa-

// [[Rcpp::export]]
IntegerMatrix CellCounts(List x, List combos) {
  
  int x_n = x.size();
  int combos_n = combos.size();
  IntegerMatrix output(x_n, combos_n);
  
  for (int i=0; i < x_n; ++i) {
    // Rcout << "Working with matrix " << i+1 << " of " << x_n << std::endl;
    LogicalMatrix mat = as<LogicalMatrix>(x[i]);
    int nrows = mat.nrow();
    for (int k=0; k < combos_n; ++k) {
      
      int num = 0;
      IntegerVector c_combo = as<IntegerVector>( combos[k] );
      int n_c = c_combo.size();
      
      IntegerVector c_combo_abs = sapply(c_combo, ::abs);
      
      for (int j=0; j < nrows; ++j) {
        
        LogicalMatrix::Row row = mat(j, _);
      
        // checking algorithm:
        // we loop through each entry in 'c_combo' => p
        // we check the element in 'row' at 'abs(c_combo[p])'
        // if 'abs(c_combo[p])' is > 0 and row[ c_combo[p] ] is true; continue
        // if 'abs(c_combo[p])' is <= 0 and row[ c_combo[p] ] is false; continue
        // else, 'false_case'
        
        #define c (c_combo[p])
        #define abs_c (c_combo_abs[p])
          
        for (int p=0; p < n_c; ++p) {
          if ((c > 0 && row[abs_c-1] <= 0) || (c < 0 && row[abs_c-1] > 0)) {
            goto false_case;
          }
        }
        
        #undef c
        #undef abs_c
        
        // if we reached here, we matched all of our conditions
        ++num;
        
        // false case allows us to skip incrementing num
        false_case: {}
        
      }
      
      // insert result into matrix
      output(i, k) = num;
      
    }
  }
  
  // set row, column names from the data if possible
  output.attr("dimnames") = List(2);
  if (!Rf_isNull(x.attr("names"))) {
    SET_VECTOR_ELT(output.attr("dimnames"), 0, Rf_duplicate(x.attr("names")));
  }
  
  if (!Rf_isNull(combos.attr("names"))) {
    SET_VECTOR_ELT(output.attr("dimnames"), 1, Rf_duplicate(combos.attr("names")));
  }
  
  return output;
  
}

/*** R
## simulate some data
K <- 6 ## number of markers
data <- replicate(10, simplify=FALSE, {
  m <- matrix( rnorm(1E4 * K, 2000, 1000 ), ncol=K )
  m[m < 2500] <- 0
  return(m)
})
names(data) <- sample(letters, 10)

combos <- list(1, 2, 3, 4, 5, 6) ## marginal cell counts
cc <- cell_counts(data, combos)
f <- function(data) {
  do.call(rbind, lapply(data, function(x) apply(x, 2, function(x) sum(x > 0))))
}
cc2 <- f(data)
identical(cc, cc2)
library(microbenchmark)
microbenchmark(
  cell_counts(data, combos),
  f(data)
)
*/
