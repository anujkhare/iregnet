/* TODO:
 * assume that sigma is fixed initially
 */

#include "iregnet.h"

// Need to take as input params
#define eps 1
#define tol 1

static  IREG_DIST get_dist_ind(Rcpp::String);

// [[Rcpp::export]]
Rcpp::NumericVector fit_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y,
                            Rcpp::String family, double alpha)
{

  std::cout<< X.ncol() << " " << X.nrow() << std::endl;
  std::cout<< y.ncol() << " " << y.nrow() << std::endl;

  for (int j=0; j < y.ncol(); ++j) {
    for (int i = 0; i < y.nrow(); ++i) {
      std::cout<<y(j, i) << " ";
    }
    std::cout<<std::endl;
  }

  std::cout << family.get_cstring() << std::endl;
  std::cout << alpha << std::endl;

  IREG_DIST dist_ind = IREG_DIST_UNKNOWN;
  dist_ind = get_dist_ind(family);
  std::cout << "family " << dist_ind << std::endl;

  /* call helper function to get the status of censoring */

  /* Initialize the parameter values using lambda_max */

  // beta should be zero at lambda_max

  // set eta = X' beta


  /* Compute the solution! */
  //return x * 2;
  return 0;
}

/*
 * Takes as input Rcpp string family name, and returns corresponding enum val
 */
static IREG_DIST get_dist_ind(Rcpp::String dist_str) {
  if (strcmp("gaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_GAUSSIAN;
  if (strcmp("logistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOGISTIC;

  return IREG_DIST_UNKNOWN;
}

/*** R

*/
