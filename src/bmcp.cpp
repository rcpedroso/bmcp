#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;






////////////////////////////////////////////////////////////////////////////////////////////////
/// runif_cpp
////////////////////////////////////////////////////////////////////////////////////////////////
double runif_cpp() {
  RNGScope scope;
  return Rf_runif(0.0,1.0);
}



////////////////////////////////////////////////////////////////////////////////////////////////
/// mu_sample
////////////////////////////////////////////////////////////////////////////////////////////////
arma::rowvec bmcp_mu_sample(arma::rowvec X, arma::rowvec U, arma::rowvec s2, double mu0, double s02) {

  int last = X.size() - 1;
  arma::rowvec mu(X.size());
  int ki = 0; int kj = 0;
  double m;
  
  RNGScope scope;

  while (1){
    while (U(kj)==1) {kj=kj+1;}
    int wr = 1;
    // int nrr = kj - ki + 1;
    arma::rowvec s2_Mr = s2.subvec(ki,kj);
    arma::rowvec X_Mr = X.subvec(ki,kj);
    double M1 = sum(s2_Mr)         +   1/(wr*s02);
    double M2 = sum(X_Mr % s2_Mr) + mu0/(wr*s02);
    m = Rf_rnorm( M2/M1 , 1/sqrt(M1) );
    for (int i = ki ; i <= kj ; i++) {mu(i) = m;}
    if (kj == last) break;
    kj = kj + 1;
    ki = kj;
  }

  // return
  return mu;
}



////////////////////////////////////////////////////////////////////////////////////////////////
/// s2_sample
////////////////////////////////////////////////////////////////////////////////////////////////
arma::rowvec bmcp_s2_sample(arma::rowvec X, arma::rowvec V, arma::rowvec mu, double a, double d) {

  int last = X.size() - 1;
  arma::rowvec s2(X.size());
  int ki = 0; int kj = 0;
  double s;
  
  RNGScope scope;

  while (1){
    while (V(kj)==1) {kj=kj+1;}
    int ns = kj - ki + 1;
    arma::rowvec mu_Ls = mu.subvec(ki,kj);
    arma::rowvec X_Ls = X.subvec(ki,kj);
    double DD = ns + d;
    double AA = sum((X_Ls-mu_Ls) % (X_Ls-mu_Ls)) + a;
    s = 1/R::rgamma( DD/2 , 2/AA );
    for (int i = ki ; i <= kj ; i++) {s2[i] = s;}
    if (kj == last) break;
    kj = kj + 1;
    ki = kj;
  }

  // return
  return s2;
}



////////////////////////////////////////////////////////////////////////////////////////////////
/// U_sample
////////////////////////////////////////////////////////////////////////////////////////////////

//////// i = 0
int bmcp_Ru1(arma::rowvec X, arma::rowvec U, int i, arma::rowvec s2, double pm, double mu0, double s02) {

  int last = X.size() - 1;

  /// i0 <- first index (change point) of cluster Mr
  int i0 = i;

  /// i1 <- last index (end point) of cluster Mr
  int k = 1;
  while (U(i+k) == 1) {
    k = k + 1;
    if (i+k == last) break;
  }
  int i1 = i + k;


  /// cluster Mr ###########
  arma::rowvec s2_Mr = s2.subvec(i0,i1);
  arma::rowvec X_Mr = X.subvec(i0,i1);
  int wr = 1;
  double Q1r = sum(s2_Mr)        +   1/(wr*s02);
  double Q2r = sum(X_Mr % s2_Mr) + mu0/(wr*s02);


  /// cluster Mri ###########
  arma::rowvec s2_Mri = s2.subvec(i0,i);
  arma::rowvec X_Mri = X.subvec(i0,i);
  int wri = 1;
  double Q1ri = sum(s2_Mri)         +   1/(wri*s02);
  double Q2ri = sum(X_Mri % s2_Mri) + mu0/(wri*s02);


  /// cluster Mir ###########
  arma::rowvec s2_Mir = s2.subvec(i+1,i1);
  arma::rowvec X_Mir = X.subvec(i+1,i1);
  int wir = 1;
  double Q1ir = sum(s2_Mir)         +   1/(wir*s02);
  double Q2ir = sum(X_Mir % s2_Mir) + mu0/(wir*s02);


  /// return
  double unif = runif_cpp();
  long double ratio = ((1-pm)/pm) /
    sqrt( exp(
        ( log(wr) + log(Q1r) - log(wri) - log(Q1ri) - log(wir) - log(Q1ir) - log(s02) ) +
          ( sum((X_Mr % X_Mr) % s2_Mr)    + (mu0*mu0)/(wr*s02)  - (Q2r*Q2r)/Q1r    ) -
          ( sum((X_Mri % X_Mri) % s2_Mri) + (mu0*mu0)/(wri*s02) - (Q2ri*Q2ri)/Q1ri ) -
          ( sum((X_Mir % X_Mir) % s2_Mir) + (mu0*mu0)/(wir*s02) - (Q2ir*Q2ir)/Q1ir )
    ) );
  return ratio >= (unif/(1-unif));
}


//////// i = 1,...,n
int bmcp_Ru(arma::rowvec X, arma::rowvec U, int i, arma::rowvec s2, double pm, double mu0, double s02) {

  int last = X.size() - 1;

  // i0 <- first index (change point) of cluster Mr
  int k = 1;
  while (U(i-k) == 1) {
    k = k + 1;
    if (i-k == -1) break;
  }
  int i0 = i - k + 1;

  // i1 <- last index (end point) of cluster Mr
  k = 1;
  while (U(i+k) == 1) {
    k = k + 1;
    if (i+k == last) break;
  }
  int i1 = i + k;


  /// cluster Mr ###########
  arma::rowvec s2_Mr = s2.subvec(i0,i1);
  arma::rowvec X_Mr = X.subvec(i0,i1);
  int wr = 1;
  double Q1r = sum(s2_Mr)        +   1/(wr*s02);
  double Q2r = sum(X_Mr % s2_Mr) + mu0/(wr*s02);


  /// cluster Mri ###########
  arma::rowvec s2_Mri = s2.subvec(i0,i);
  arma::rowvec X_Mri = X.subvec(i0,i);
  int wri = 1;
  double Q1ri = sum(s2_Mri)         +   1/(wri*s02);
  double Q2ri = sum(X_Mri % s2_Mri) + mu0/(wri*s02);


  /// cluster Mir ###########
  arma::rowvec s2_Mir = s2.subvec(i+1,i1);
  arma::rowvec X_Mir = X.subvec(i+1,i1);
  int wir = 1;
  double Q1ir = sum(s2_Mir)         +   1/(wir*s02);
  double Q2ir = sum(X_Mir % s2_Mir) + mu0/(wir*s02);


  /// return
  double unif = runif_cpp();
  long double ratio = ((1-pm)/pm) /
    sqrt( exp(
        ( log(wr) + log(Q1r) - log(wri) - log(Q1ri) - log(wir) - log(Q1ir) - log(s02) ) +
          ( sum((X_Mr % X_Mr) % s2_Mr)    + (mu0*mu0)/(wr*s02)  - (Q2r*Q2r)/Q1r    ) -
          ( sum((X_Mri % X_Mri) % s2_Mri) + (mu0*mu0)/(wri*s02) - (Q2ri*Q2ri)/Q1ri ) -
          ( sum((X_Mir % X_Mir) % s2_Mir) + (mu0*mu0)/(wir*s02) - (Q2ir*Q2ir)/Q1ir )
    ) );
  return ratio >= (unif/(1-unif));
}



////////////////////////////////////////////////////////////////////////////////////////////////
/// V_sample
////////////////////////////////////////////////////////////////////////////////////////////////

//////// i = 0
int bmcp_Rv1(arma::rowvec X, arma::rowvec V, int i, arma::rowvec mu, double ps, double a, double d) {

  int last = V.size() - 1;

  /// i0 <- first index (change point) of cluster Ls
  int i0 = i;

  /// i1 <- last index (end point) of cluster Ls
  int k = 1;
  while (V(i+k) == 1) {
    k = k + 1;
    if (i+k == last) break;
  }
  int i1 = i + k;


  /// cluster Ls ###########
  int ns = i1 - i0 + 1;
  arma::rowvec mu_Ls = mu.subvec(i0,i1);
  arma::rowvec X_Ls = X.subvec(i0,i1);
  double As = sum((X_Ls-mu_Ls) % (X_Ls-mu_Ls)) + a;
  double Ds = ns + d;


  /// cluster Lsi ###########
  int nsi = i - i0 + 1;
  arma::rowvec mu_Lsi = mu.subvec(i0,i);
  arma::rowvec X_Lsi = X.subvec(i0,i);
  double Asi = sum((X_Lsi-mu_Lsi) % (X_Lsi-mu_Lsi)) + a;
  double Dsi = nsi + d;


  /// cluster Lis ###########
  int nis = i1 - i;
  arma::rowvec mu_Lis = mu.subvec(i+1,i1);
  arma::rowvec X_Lis = X.subvec(i+1,i1);
  double Ais = sum((X_Lis-mu_Lis) % (X_Lis-mu_Lis)) + a;
  double Dis = nis + d;


  /// return
  double unif = runif_cpp();
  double ratio = (1-ps)/ps * tgamma(d/2)/pow(a/2,d/2) *
    // exp( lgamma(Ds/2) - lgamma(Dsi/2) - lgamma(Dis/2) ) *
    // exp( (Dsi/2)*log(Asi/2) + (Dis/2)*log(Ais/2) - (Ds/2)*log(As/2) );
    exp( lgamma(Ds/2) - lgamma(Dsi/2) - lgamma(Dis/2) +
    (Dsi/2)*log(Asi/2) + (Dis/2)*log(Ais/2) - (Ds/2)*log(As/2) );
  return ratio >= (unif/(1-unif));
}


//////// i = 1,...,n
int bmcp_Rv(arma::rowvec X, arma::rowvec V, int i, arma::rowvec mu, double ps, double a, double d) {

  int last = V.size() - 1;

  // i0 <- first index (change point) of cluster Ls
  int k = 1;
  while (V(i-k) == 1) {
    k = k + 1;
    if (i-k == -1) break;
  }
  int i0 = i - k + 1;

  // i1 <- last index (end point) of cluster Ls
  k = 1;
  while (V(i+k) == 1) {
    k = k + 1;
    if (i+k == last) break;
  }
  int i1 = i + k;


  /// cluster Ls ###########
  int ns = i1 - i0 + 1;
  arma::rowvec mu_Ls = mu.subvec(i0,i1);
  arma::rowvec X_Ls = X.subvec(i0,i1);
  double As = sum((X_Ls-mu_Ls) % (X_Ls-mu_Ls)) + a;
  double Ds = ns + d;


  /// cluster Lsi ###########
  int nsi = i - i0 + 1;
  arma::rowvec mu_Lsi = mu.subvec(i0,i);
  arma::rowvec X_Lsi = X.subvec(i0,i);
  double Asi = sum((X_Lsi-mu_Lsi) % (X_Lsi-mu_Lsi)) + a;
  double Dsi = nsi + d;


  /// cluster Lis ###########
  int nis = i1 - i;
  arma::rowvec mu_Lis = mu.subvec(i+1,i1);
  arma::rowvec X_Lis = X.subvec(i+1,i1);
  double Ais = sum((X_Lis-mu_Lis) % (X_Lis-mu_Lis)) + a;
  double Dis = nis + d;

  /// return
  double unif = runif_cpp();
  double ratio = (1-ps)/ps * tgamma(d/2)/pow(a/2,d/2) *
    // exp( lgamma(Ds/2) - lgamma(Dsi/2) - lgamma(Dis/2) ) *
    // exp( (Dsi/2)*log(Asi/2) + (Dis/2)*log(Ais/2) - (Ds/2)*log(As/2) );
    exp( lgamma(Ds/2) - lgamma(Dsi/2) - lgamma(Dis/2) +
    (Dsi/2)*log(Asi/2) + (Dis/2)*log(Ais/2) - (Ds/2)*log(As/2) );
  return ratio >= (unif/(1-unif));
}



////////////////////////////////////////////////////////////////////////////////////////////////
/// bmcp
////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List bmcp(int burn, int N,               // sample configuration
                arma::rowvec X,                // observations
                double alpha1, double beta1,   // p1 prior
                double alpha2, double beta2,   // p2 prior
                double a, double d,            // s2 prior
                double mu0, double s02) {      // mu prior
  
  int n = X.size();
  arma::Col<int> b1(N), b2(N);
  arma::colvec p1(N), p2(N);
  arma::mat s2 = arma::zeros(N,n);
  arma::mat mu = arma::zeros(N,n);
  arma::mat u = arma::zeros(N,n);
  arma::mat v = arma::zeros(N,n);
  
  RNGScope scope;
  
  // Initial values
  p1(0) = alpha1/(alpha1 + beta1);
  p2(0) = alpha2/(alpha2 + beta2);
  b1(0) = n - sum(u.row(0));
  b2(0) = n - sum(v.row(0));
  s2.row(0) = (X - mean(X)) % (X - mean(X));
  mu.row(0) = X;
  

  for(int s = 1; s < N; s++) {

    // p sample
    p1(s) = Rf_rbeta( alpha1+b1(s-1)-1 , n+beta1-b1(s-1) );
    p2(s) = Rf_rbeta( alpha2+b2(s-1)-1 , n+beta2-b2(s-1) );
    // Rprintf("the value of p1 : %f \n", p1(s));
    // Rprintf("the value of p2 : %f \n", p2(s));
    

    // U sample
    u.row(s) = u.row(s-1);
    u(s,0) = bmcp_Ru1(X , u.row(s) , 0 , 1/s2.row(s-1) , p1(s), mu0 , s02);
    // Rprintf("the value of u(s,0) : %i \n", u(s,0));
    for(int pos = 1; pos < n-1; pos++) {
      u(s,pos) = bmcp_Ru(X , u.row(s) , pos , 1/s2.row(s-1) , p1(s), mu0 , s02);
    }
    // Rprintf("the value of u(s,1) : %i \n", u(s,1));
    
    
    // mu sample
    mu.row(s) = bmcp_mu_sample(X , u.row(s) , 1/s2.row(s-1) , mu0 , s02);
    // Rprintf("the value of mu(s,0) : %f \n", mu(s,0));


    // V sample
    v.row(s) = v.row(s-1);
    v(s,0) = bmcp_Rv1(X , v.row(s) , 0 , mu.row(s) , p2(s) , a , d);
    // Rprintf("the value of v(s,0) : %i \n", u(s,0));
    for(int pos = 1; pos < n-1; pos++) {
      v(s,pos) = bmcp_Rv(X , v.row(s) , pos , mu.row(s) , p2(s) , a , d);
    }
    // Rprintf("the value of v(s,1) : %i \n", v(s,1));


    // s2 sample
    s2.row(s) = bmcp_s2_sample(X , v.row(s) , mu.row(s) , a , d);
    // Rprintf("the value of s2(s,0) : %f \n", s2(s,0));


    // Number of endpoints/blocks
    b1(s) = n - sum(u.row(s));
    b2(s) = n - sum(v.row(s));
    // Rprintf("the value of b1(s) : %f \n", b1(s));
    // Rprintf("the value of b2(s) : %f \n", b2(s));
    // Rprintf("the value of s : %i \n", s);

  }

  List out;
  out["p1"] = p1.subvec( burn, N-1 );
  out["p2"] = p2.subvec( burn, N-1 );
  out["b1"] = b1.subvec( burn, N-1 );
  out["b2"] = b2.subvec( burn, N-1 );
  out["mu"] = mu.submat( burn, 0, N-1, n-1 );
  out["s2"] = s2.submat( burn, 0, N-1, n-1 );;
  out["u"] = u.submat( burn, 0, N-1, n-1 );;
  out["v"] = v.submat( burn, 0, N-1, n-1 );;
  out["alpha1"] = alpha1;
  out["beta1"] = beta1;
  out["alpha2"] = alpha2;
  out["beta2"] = beta2;
  out["mu0"] = mu0;
  out["s02"] = s02;
  out["a"] = a;
  out["d"] = d;

  return out;
}







