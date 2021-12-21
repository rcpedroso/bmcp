#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix u_to_index(NumericMatrix U) {
  
  int nc = U.ncol();
  int nr = U.nrow();
  
  NumericMatrix eps(nr,nc+1);

  for (int k = 0; k < nr; k++) {
  
    int i = -1;
  
    for (int j = 0; j < nc; j++) {
      if (U(k,j)==0) {
        i = i + 1;
        eps(k,i) = j + 1;
      }
      eps(k,nc) = i + 1;
    }
  }
  return(eps);
}

 
  
