#include <RcppArmadillo.h>
#include "prox.h"

using namespace Rcpp;
using namespace arma;

//' Sorted L1 Norm Prox
//'
//' @param x coefficients
//' @param x regularization
//'
//' @return Prox
//'
//' @export
// [[Rcpp::export]]
arma::vec sorted_l1_prox(arma::vec x, arma::vec lambda)
{
  const uword p = x.n_elem;

  Slope prox(p);

  return prox(x, lambda);
}

