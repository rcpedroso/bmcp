#ifndef BOUNDS_H
#define BOUNDS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


/* last position of Sj (for i=0) */
int sj1(arma::rowvec U, int last) {
  
  /* i1 <- last index (end point) of cluster Sj */
  int i1 = 1;
  while (U[i1] == 1) {
    i1 += 1;
    if (i1 == last) break;
  }

  // return
  return i1;
}



/* first and last position of Sj (for i>0) */
arma::Col<int> sj01(arma::rowvec U, int last, int i) {
  
  /* i0 <- first index (change point) of cluster Sj */
  int k = 1;
  while (U[i-k] == 1) {
    k = k + 1;
    if (i-k == -1) break;
  }
  int i0 = i - k + 1;

  /* i1 <- last index (end point) of cluster Sj */
  k = 1;
  while (U[i+k] == 1) {
    k = k + 1;
    if (i+k == last) break;
  }
  int i1 = i + k;

  // return
  return {i0,i1};
}

  
#endif