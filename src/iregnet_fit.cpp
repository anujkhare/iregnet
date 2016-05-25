/* TODO:
 * assume that sigma is fixed initially
 */
// used for all counts
#define ull unsigned long long

#include "iregnet.h"

static IREG_DIST get_ireg_dist(Rcpp::String);
static inline void get_censoring_types (Rcpp::NumericMatrix y, IREG_CENSORING *censoring_type);

/* fit_cpp: Fit a censored data distribution with elastic net reg.
 *
 * Inputs:
 *
 * Outputs:
 *      beta:     the final coef vector
 *      score:    the final score vector
 *      loglik:   the final log-likelihood
 *      n_iters:  number of iterations consumed
 *      error_status: 0 = no errors, -1 = input error, -2 = did not converge
 *
 * Work arrays created:
 *      ?
 */
// [[Rcpp::export]]
Rcpp::List fit_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y,
                            Rcpp::String family,   double alpha,
                            bool estimate_scale = false, double tol_chol = 0.1,
                            double maxiter = 100,  double tol_convergence = 0.1,
                            int flag_debug = 0)
{

  /* Uselesss print stuff for now */
  if (flag_debug == IREG_DEBUG_INPUT) {
    // Rcpp has implemented << for NumericVectors and matrices too!
    std::cout << "X:\n" << X << std::endl;
    std::cout << "y:\n" << y << std::endl;
    std::cout << "family:" << family.get_cstring() << std::endl;
    std::cout << "alpha:" << alpha << std::endl;
  }


  /* Initialise some helper variables */
  ull n_obs, n_vars, n_cols_y, n_params;
  IREG_DIST dist;

  n_obs  = X.nrow();
  n_cols_y = y.ncol();
  dist = get_ireg_dist(family);

  n_vars = X.ncol();  // n_vars is the number of variables corresponding to the coeffs of X
  n_params = n_vars + int(estimate_scale); // n_params is the number of parameters
                                           // to optimize (includes sigma as well)

  if (flag_debug == IREG_DEBUG_N) {
    std::cout << "n_vars: " << n_vars << ", n_params: " << n_params << "\n";
    std::cout << "n_obs: " << n_obs << "\n";
  }


  /* Validate all the arguments again */
  if ((alpha > 1 || alpha < 0) || (IREG_DIST_UNKNOWN == dist) ||
      (y.nrow() != n_obs)) {
    //error_status = -1;      // TODO: UNUSED!
    return Rcpp::List::create(Rcpp::Named("error_status") = -1);
  }


  /* Create output variables */
  Rcpp::NumericVector out_beta(n_params);
  Rcpp::NumericVector out_score(n_params);
  double *beta  = &out_beta[0];     // pointers for the C++ routines
  double *score = &out_score[0];
  double loglik;
  ull    n_iters;


  /* get censoring types of the observations */
  IREG_CENSORING censoring_type[n_obs];
  get_censoring_types(y, censoring_type); // NANs denote censored observations in y

  if (flag_debug == IREG_DEBUG_CENSORING) {
    std::cout << "Censoring types:\n";
    for (int i = 0; i < n_obs; ++i)
      std::cout << censoring_type[i] << std::endl;
  }

  /* Initialize the parameter values using lambda_max */
  // beta should be zero at lambda_max
  // set eta = X' beta
  /* Compute the solution! */

  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta,
                            Rcpp::Named("score")        = out_score,
                            Rcpp::Named("loglik")       = loglik,
                            Rcpp::Named("error_status") = 0,
                            Rcpp::Named("n_iters")      = n_iters);
}


static inline void get_censoring_types (Rcpp::NumericMatrix y, IREG_CENSORING *censoring_type)
{
  for (ull i = 0; i < y.nrow(); ++i) {
    // std::cout << y(i, 0) << " " << y(i, 1) << "\n";
    if (y(i, 0) == Rcpp::NA) {

      if (y(i, 1) == Rcpp::NA)
        censoring_type[i] = IREG_CENSOR_INVALID;    // invalid data
      else
        censoring_type[i] = IREG_CENSOR_LEFT;       // left censoring

    } else {
      if (y(i, 1) == Rcpp::NA)                      // right censoring
        censoring_type[i] = IREG_CENSOR_RIGHT;

      else if (y(i, 0) == y(i, 1))
          censoring_type[i] = IREG_CENSOR_NONE;     // no censoring
        else
          censoring_type[i] = IREG_CENSOR_INTERVAL; // interval censoring
    }
  } // end for
}


/*
 * Takes as input Rcpp string family name, and returns corresponding enum val
 */
static IREG_DIST get_ireg_dist (Rcpp::String dist_str) {
  if (strcmp("gaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_GAUSSIAN;
  if (strcmp("logistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOGISTIC;

  return IREG_DIST_UNKNOWN;
}

/*** R

*/
