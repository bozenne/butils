// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
void hello0() {
  Rcpp::Rcout << "Hello!\n";
}

// [[Rcpp::export]]
arma::mat hello1(arma::mat A) {
  Rcpp::Rcout << "Hello1!\n";
  return A;
}

// [[Rcpp::export]]
arma::mat hello2(arma::mat A, SEXP B) {
  Rcpp::Rcout << "Hello2!\n";
  return A;
}