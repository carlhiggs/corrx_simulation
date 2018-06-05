#include <Rcpp.h>
#include <numeric>

// [[Rcpp::export]]
Rcpp::NumericVector rnorm_cpp(int n, double m, double s){
	Rcpp::RNGScope scope;
    Rcpp::NumericVector y = Rcpp::rnorm(n,m,s);
    return(y);
    }

// [[Rcpp::export]]
Rcpp::NumericVector rchisq_cpp(int n, double df, double ncp){
	Rcpp::RNGScope scope;
    Rcpp::NumericVector y = Rcpp::rnchisq(n,df,ncp);
    return(y);
    }

