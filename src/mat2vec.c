#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>

void set_rownames(SEXP x);
void set_names(SEXP x);

SEXP mat2vec(SEXP x) {
  char type = TYPEOF(x);
  if (!isMatrix(x))
    error("'x' must be a matrix");
  int nRow = INTEGER(getAttrib(x, R_DimSymbol))[0];
  int nCol = INTEGER(getAttrib(x, R_DimSymbol))[1];
  SEXP output = PROTECT( allocVector(VECSXP, nCol) );
  
  #define HANDLE_CASE(ACCESSOR, CTYPE) { \
  for (int i=0; i < nCol; ++i) { \
    SET_VECTOR_ELT(output, i, allocVector(type, nRow)); \
    CTYPE* output_ptr = ACCESSOR(VECTOR_ELT(output, i)); \
    CTYPE* x_ptr = ACCESSOR(x); \
    for (int j=0; j < nRow; ++j) { \
      output_ptr[j] = x_ptr[nRow*i + j]; \
    } \
  } \
  break; \
} \

#define HANDLE_CASE_MEMCPY(ACCESSOR, CTYPE) { \
for (int i=0; i < nCol; ++i) { \
  char sz = sizeof(CTYPE); \
  SET_VECTOR_ELT(output, i, allocVector(type, nRow)); \
  SEXP elt = VECTOR_ELT(output, i); \
  memcpy( \
    (char*) DATAPTR(elt), \
    (char*) DATAPTR(x) + (i*nRow*sz), \
    nRow*sz \
  ); \
} \
break; \
} \

switch (type) {
  case INTSXP: HANDLE_CASE_MEMCPY(INTEGER, int);
  case REALSXP: HANDLE_CASE_MEMCPY(REAL, double);
  case LGLSXP: HANDLE_CASE_MEMCPY(INTEGER, int);
  case STRSXP: HANDLE_CASE(STRING_PTR, SEXP);
  default: error("Unhandled SEXP type '%s'", type2char(type));
}

#undef HANDLE_CASE
#undef HANDLE_CASE_MEMCPY

setAttrib(output, R_NamesSymbol, VECTOR_ELT(getAttrib(x, R_DimNamesSymbol), 1));

UNPROTECT(1);
return output;

}

#undef USE_RINTERNALS
