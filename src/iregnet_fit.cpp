/* TODO:
 * assume that sigma is fixed initially
 */

#include "iregnet.h"

#define EPSILON_LAMBDA 0.0001
#define M_LAMBDA  1

static IREG_DIST get_ireg_dist(Rcpp::String);
static inline void get_censoring_types (Rcpp::NumericMatrix y, IREG_CENSORING *censoring_type);
static inline double soft_threshold(double x, double lambda);

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
                   bool estimate_scale = false, double scale = 0,
                   double tol_chol = 0.1,
                   double maxiter = 100,  double tol_convergence = 0.1,
                   int flag_debug = 0)
{

  /* Uselesss print stuff for now */
  if (flag_debug == IREG_DEBUG_INPUT) {       // 1
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
                                           // to optimize (includes scale as well)

  if (flag_debug == IREG_DEBUG_N) {         // 2
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
  //Rcpp::NumericVector out_score(n_params);
  //Rcpp::NumericMatrix info_mat(n_params, n_params);
  //Rcpp::NumericMatrix var_mat(n_params, n_params);
  double *beta  = REAL(out_beta);     // pointers for the C++ routines
  //double *score_ptr = REAL(out_score);
  double loglik;
  ull    n_iters;

  // Temporary variables: not returned // TODO: Maybe alloc them together?
  double *eta = new double [n_obs];         // vector of linear predictors = X' beta
  // eta = 0 for the initial lambda_max, and in each iteration of coordinate descent,
  // eta is updated along with beta in place

  //double *mu = new double [n_obs];          // grad of LL wrt eta
  double *w  = new double [n_obs];          // diagonal of the hessian of LL wrt eta
                                            // these are the weights of the IRLS linear reg
  double *z = new double [n_obs];           // z_i = eta_i - mu_i / w_i
  double lambda_max, lambda_min, n_lambda;  // for the pathwise solution

  //double *arr = new double [n_params * n_params + n_params * n_params];
  //double **info_mat_ptr = dmatrix(arr, n_params, n_params);
  //double **var_mat_ptr = dmatrix(arr + n_params * n_params, n_params, n_params);
  //double **x_ptr = dmatrix(REAL(X), n_vars, n_obs);

  /* get censoring types of the observations */
  IREG_CENSORING censoring_type[n_obs];
  get_censoring_types(y, censoring_type); // NANs denote censored observations in y

  if (flag_debug == IREG_DEBUG_CENSORING) {         // 3
    Rcpp::NumericVector rvec (censoring_type, censoring_type + n_obs);
    return Rcpp::List::create(Rcpp::Named("error_status")   = 0,
                              Rcpp::Named("censoring_types") = rvec);
  }


  /* Apply the required transformation to the time scale (t = log(y)) */
  // std::cout << y << std::endl;
  for (ull i = 0; i < 2 * n_obs; ++i) {     // Rcpp::NumericMatrix is unrolled col major
    //std::cout << y[i] << std::endl;
    y[i] = log(y[i]);
  }


  /* Initalize the pathwise solution
   * We will always start with lambda_max, at which beta = 0, eta = 0.
   * If some value of lambda is supplied, we will stop iterations at that value.
   */

  // TODO: FIGURE OUT WHAT HAPPENS TO THE SIGMA TERMS HERE
  // set beta = 0 and eta = 0
  for (ull i = 0; i < n_params; ++i) {
    beta[i] = eta[i] = 0;
  }

  /////////////////////////////////////////
  // Calculate w and z right here!
  compute_grad_response(w, z, REAL(y), REAL(y) + n_obs, eta, scale,
                        censoring_type, n_obs, dist);

  // Calculate lambda_max
    // TODO: try to optimize by reversing j and i loops and using an extra array
  lambda_max = -1;
  for (ull j = 0; j < n_params; ++j) {
    double temp = 0;

    for (ull i = 0; i < n_obs; ++i) {
      temp += (w[i] * X(i, j) * z[i]);
    }
    //temp = temp / (alpha);

    lambda_max = (lambda_max > temp)? lambda_max: temp;
  }

  lambda_min = EPSILON_LAMBDA * lambda_max;


  /* Iterate over grid of lambda values */
  double lambda, lambda_ratio;
  bool flag_beta_converged = 0;
  //double w_x_z = new double [n_obs];
  //double *w_x_eta = new double [n_obs];
  double w_x_z, w_x_eta, sol_num, sol_denom;

  lambda_ratio = (lambda_max / lambda_min);

  for (int m = M_LAMBDA; m >= 0; --m) {
    lambda = std::pow(lambda_ratio, (1.0 * m) / M_LAMBDA);


    // calculate the intermediate terms in the soft threshold expr
    // whether or not to store this depends on how many times the beta iterations run
    //w_x_z = 0;
    //for (ull i = 0; i < n_obs; ++i) {
    //  w_x_z += (w[i] * X(i, k)
    //}

    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    flag_beta_converged = 0;              // =1 if beta converges
    do {                                  // until Convergence of beta

      /* iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_params; ++k) {

        sol_num = sol_denom = 0;          // TODO: You should optimize this so that we don't calculate the whole thing everytime
        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
          sol_num += (w[i] * X(i, k) * (z[i] - eta[i]));
          sol_denom += (w[i] * X(i, k) * X(i, k));
        }

        beta[k] = soft_threshold(sol_num / (sol_denom + lambda * (1 - alpha)), lambda * alpha);

        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
        }

      }   // end for: beta_k solution

      flag_beta_converged = 1;    // TODO: CHECK FOR CONVERGENCE OF BETA! (or max iters for beta)
    } while (flag_beta_converged != 1);

    /* beta and eta will already contain their updated values since we calculate them in place */
    // calculate w and z again (beta & hence eta would have changed)
    compute_grad_response(w, z, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                          censoring_type, n_obs, dist);

  } // end for: lambda


  /* Free the temporary variables */
  delete [] eta;
  delete [] w;
  delete [] z;

  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta,
                            //Rcpp::Named("score")        = out_score,
                            Rcpp::Named("loglik")       = loglik,
                            Rcpp::Named("error_status") = 0,
                            Rcpp::Named("n_iters")      = n_iters);
}

  /* Convert the input to normalized form e_i = (log(y_i) - x_i' beta) / sigma */

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
static IREG_DIST get_ireg_dist (Rcpp::String dist_str)
{
  if (strcmp("gaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_GAUSSIAN;
  if (strcmp("logistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOGISTIC;

  return IREG_DIST_UNKNOWN;
}

static inline double soft_threshold(double x, double lambda)
{
  double temp = abs(x) - lambda;
  return (temp > 0)? ((x > 0)? 1: -1) * temp: 0;
}
