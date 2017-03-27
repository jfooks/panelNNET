#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]

cppFunction('arma::mat logistic(arma::mat v, double s = 1) {
  return 1/(1+exp(-v*s));
}
', depends="RcppArmadillo")



