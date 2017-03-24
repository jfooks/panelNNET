#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector logistic(NumericVector v, double s = 1) {
  return 1/(1+exp(-v*s));
}

