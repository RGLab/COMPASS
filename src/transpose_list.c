#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>

SEXP transpose_list(SEXP x) {
  
  int n = length(x);
  int N = length(VECTOR_ELT(x, 0));
  int type = TYPEOF( VECTOR_ELT(x, 0) );
  
  SEXP tmp;
  SEXP output = PROTECT( allocVector(VECSXP, N) );
  
  #define HANDLE_CASE(RTYPE, CTYPE, ACCESSOR) { \
    tmp = PROTECT( allocVector(RTYPE, n) ); \
    CTYPE* ptr = ACCESSOR(tmp); \
    for (int j=0; j < N; ++j ) { \
      for (int i=0; i < n; ++i) { \
        ptr[i] = ACCESSOR( VECTOR_ELT(x, i) )[j]; \
      } \
      SET_VECTOR_ELT(output, j, duplicate(tmp)); \
    } \
    UNPROTECT(1); \
    break; \
  } \
    
  switch (type) {
  case LGLSXP: HANDLE_CASE(LGLSXP, int, LOGICAL);
  case INTSXP: HANDLE_CASE(INTSXP, int, INTEGER);
  case REALSXP: HANDLE_CASE(REALSXP, double, REAL);
  case STRSXP: HANDLE_CASE(STRSXP, SEXP, STRING_PTR);
  }
  
  UNPROTECT(1);
  return output;
  
}

#undef USE_RINTERNALS
