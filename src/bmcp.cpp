#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
#include "bounds.h"



////* mu_sample *////
arma::rowvec bmcp_mu_sample(arma::rowvec X, arma::rowvec U, arma::rowvec s2, double mu0, double s02) {

	int last = X.size() - 1;
	arma::rowvec mu(X.size());
	int ki = 0; int kj = 0;
	double m;

	RNGScope scope;

	while (1){
		while (U(kj)==1) {kj=kj+1;}
		arma::rowvec s2j = s2.subvec(ki,kj);
		arma::rowvec Xj = X.subvec(ki,kj);
		double M1 = sum(s2j)      +   1/s02;
		double M2 = sum(Xj % s2j) + mu0/s02;
		m = Rf_rnorm( M2/M1 , 1/sqrt(M1) );
		for (int i = ki ; i <= kj ; i++) {mu[i] = m;}
		if (kj == last) break;
		kj = kj + 1;
		ki = kj;
	}

	return mu;
}



////* s2_sample *////
arma::rowvec bmcp_s2_sample(arma::rowvec X, arma::rowvec U, arma::rowvec mu, double a, double d) {

	int last = X.size() - 1;
	arma::rowvec s2(X.size());
	int ki = 0; int kj = 0;
	double s;

	RNGScope scope;

	while (1){
		while (U(kj)==1) {kj=kj+1;}
		int nj = kj - ki + 1;
		arma::rowvec muj = mu.subvec(ki,kj);
		arma::rowvec Xj = X.subvec(ki,kj);
		double D = nj + d;
		double A = sum((Xj-muj) % (Xj-muj)) + a;
		s = 1/R::rgamma( D/2 , 2/A );
		for (int i = ki ; i <= kj ; i++) {s2[i] = s;}
		if (kj == last) break;
		kj = kj + 1;
		ki = kj;
	}

	return s2;
}



////* U1 sample *////
int bmcp_U1(arma::rowvec X, int i0, int i, int i1, arma::rowvec s2, double p1, double mu0, double s02) {

	/* cluster Sj */
	arma::rowvec s2j = s2.subvec(i0,i1);
	arma::rowvec Xj = X.subvec(i0,i1);
	double Q1j = sum(s2j)      +   1/s02;
	double Q2j = sum(Xj % s2j) + mu0/s02;

	/* cluster Sj(1) */
	arma::rowvec s2j1 = s2.subvec(i0,i);
	arma::rowvec Xj1 = X.subvec(i0,i);
	double Q1j1 = sum(s2j1)       +   1/s02;
	double Q2j1 = sum(Xj1 % s2j1) + mu0/s02;

	/* cluster Sj(2) */
	arma::rowvec s2j2 = s2.subvec(i+1,i1);
	arma::rowvec Xj2 = X.subvec(i+1,i1);
	double Q1j2 = sum(s2j2)       +   1/s02;
	double Q2j2 = sum(Xj2 % s2j2) + mu0/s02;

	/* sample */
	RNGScope scope;
	double unif = Rf_runif(0.0,1.0);
	long double ratio = ((1-p1)/p1) /
		sqrt( exp(
			( log(Q1j) - log(Q1j1) - log(Q1j2) - log(s02) - (mu0*mu0)/s02 ) +
			( sum((Xj  % Xj)  % s2j ) - (Q2j *Q2j )/Q1j  ) -
			( sum((Xj1 % Xj1) % s2j1) - (Q2j1*Q2j1)/Q1j1 ) -
			( sum((Xj2 % Xj2) % s2j2) - (Q2j2*Q2j2)/Q1j2 ) ) );
	return ratio >= (unif/(1-unif));
}



////* U2 sample *////
int bmcp_U2(arma::rowvec X, int i0, int i, int i1, arma::rowvec mu, double p2, double a, double d) {

	/* cluster Sj */
	int nj = i1 - i0 + 1;
	arma::rowvec muj = mu.subvec(i0,i1);
	arma::rowvec Xj = X.subvec(i0,i1);
	double Aj = sum((Xj-muj) % (Xj-muj)) + a;
	double Dj = nj + d;

	/* cluster Sj(1) */
	int nj1 = i - i0 + 1;
	arma::rowvec muj1 = mu.subvec(i0,i);
	arma::rowvec Xj1 = X.subvec(i0,i);
	double Aj1 = sum((Xj1-muj1) % (Xj1-muj1)) + a;
	double Dj1 = nj1 + d;

	/* cluster Sj(2) */
	int nj2 = i1 - i;
	arma::rowvec muj2 = mu.subvec(i+1,i1);
	arma::rowvec Xj2 = X.subvec(i+1,i1);
	double Aj2 = sum((Xj2-muj2) % (Xj2-muj2)) + a;
	double Dj2 = nj2 + d;

	/* sample */
	RNGScope scope;
	double unif = Rf_runif(0.0,1.0);
	double ratio = (1-p2)/p2 * tgamma(d/2)/pow(a/2,d/2) *
		exp( lgamma(Dj/2) - lgamma(Dj1/2) - lgamma(Dj2/2) +
			 (Dj1/2)*log(Aj1/2) + (Dj2/2)*log(Aj2/2) - (Dj/2)*log(Aj/2) );
	return ratio >= (unif/(1-unif));
}



////* bmcp *////
// [[Rcpp::export]]
Rcpp::List bmcp(arma::rowvec X,
                double alpha1, double beta1, double alpha2, double beta2,
                double a, double d, double mu0, double s02,
				        const int burn=1000, const int ns=1000, const int thin=5) {
  
	int n = X.size();
	int last = n - 1;
	arma::Col<int> b1(ns), b2(ns);
	arma::colvec p1(ns), p2(ns);
	arma::mat s2 = arma::zeros(ns,n);
	arma::mat mu = arma::zeros(ns,n);
	arma::mat u1 = arma::zeros(ns,n);
	arma::mat u2 = arma::zeros(ns,n);
	arma::Col<int> indexes(2);
	int i1;

	RNGScope scope;

	/* Initial values */
	p1(0) = alpha1/(alpha1 + beta1);
	p2(0) = alpha2/(alpha2 + beta2);
	b1(0) = n - sum(u1.row(0));
	b2(0) = n - sum(u2.row(0));
	s2.row(0) = (X - mean(X)) % (X - mean(X));
	mu.row(0) = X;
  

  /* mcmc */
	int iter = 0;
	int s = 0;
	while (s < ns) {
		
		iter += 1;

		/* p sample */
		p1(s) = Rf_rbeta( alpha1+b1(s)-1 , n+beta1-b1(s) );
		p2(s) = Rf_rbeta( alpha2+b2(s)-1 , n+beta2-b2(s) );
    

		/* U1 sample */
		// i = 0
		i1 = sj1(u1.row(s), last);
		u1(s,0) = bmcp_U1(X, 0, 0, i1, 1/s2.row(s), p1(s), mu0, s02);
		// i > 0
		for(int i = 1; i < n-1; i++) {
			indexes = sj01(u1.row(s), last, i);
			u1(s,i) = bmcp_U1(X, indexes(0), i, indexes(1), 1/s2.row(s), p1(s), mu0, s02);
		}
    
    
		/* mu sample */
		mu.row(s) = bmcp_mu_sample(X, u1.row(s), 1/s2.row(s), mu0, s02);


		/* U2 sample */
		// i = 0
		i1 = sj1(u2.row(s), last);
		u2(s,0) = bmcp_U2(X, 0, 0, i1, mu.row(s), p2(s), a, d);
		// i > 0
		for(int i = 1; i < n-1; i++) {
			indexes = sj01(u2.row(s), last, i);
			u2(s,i) = bmcp_U2(X, indexes(0), i, indexes(1), mu.row(s), p2(s), a, d);
		}


		/* s2 sample */
		s2.row(s) = bmcp_s2_sample(X, u2.row(s), mu.row(s), a, d);


		/* number of blocks */
		b1(s) = n - sum(u1.row(s));
		b2(s) = n - sum(u2.row(s));
		

		/* save sample */
		if (iter >= burn) {
		  if ((iter-burn) % thin == 0) {
  			s += 1;
  			if (s < ns) {
				
				p1(s) = p1(s-1);
				u1.row(s) = u1.row(s-1);
  				mu.row(s) = mu.row(s-1);
				b1(s) = b1(s-1);
				
				p2(s) = p2(s-1);
				u2.row(s) = u2.row(s-1);
  				s2.row(s) = s2.row(s-1);
				b2(s) = b2(s-1);
				
  			}
		  }
		}

	}

	List out;
  
	// data
	out["X"] = X;
	  
	// parameter chains
	out["p1"] = p1;
	out["p2"] = p2;
	out["b1"] = b1;
	out["b2"] = b2;
	out["mu"] = mu;
	out["s2"] = s2;
	out["u1"] = u1;
	out["u2"] = u2;
	
	// hyperparameters
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







