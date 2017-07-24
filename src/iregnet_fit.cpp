/* for easy compilation in emacs -*- compile-command: "R CMD INSTALL .." -*- */

/*
 * iregnet_fit.cpp
 * Author: Anuj Khare <khareanuj18@gmail.com>
 * We use 2 spaces per tab, and expand the tabs
 */
#include "iregnet.h"

#define BIG 1e35

static inline double
get_y_means (mat &y, IREG_CENSORING *status, double *ym);
static inline double soft_threshold (double x, double lambda);
static double get_init_var (double *ym, IREG_CENSORING *status, ull n, IREG_DIST dist);
static void standardize_x (mat &X,
                           double *mean_x, double *std_x,
                           bool intercept);
static void
standardize_y (mat &y, double *ym, double &mean_y);
static inline double compute_lambda_max(mat X, rowvec *w, rowvec *z,
                                        rowvec *eta, bool intercept, double &alpha,
                                        ull n_vars, ull n_obs, bool debug);

static int *
sort_X_and_y(mat &sorted_X, mat &sorted_y, IREG_CENSORING *status, mat X, mat y, int function_type);

double
identity (double y)
{
  return y;
}

static double
max(double a, double b)
{
  if(a > b)
    return a;
  return b;
}

double (*target_compute_grad_response)(rowvec *w, rowvec *z, double *scale_update, const rowvec *y_l, const rowvec *y_r,
                                       const rowvec &eta, const double scale, const IREG_CENSORING *censoring_type,
                                       const ull n_obs, IREG_DIST dist, double *mu, bool debug, const bool estimate_scale,const rowvec &y_eta,
                                       const rowvec &y_eta_square,const int *separator, rowvec *tempvar);

/* fit_cpp: Fit a censored data distribution with elastic net reg.
 * Outputs:
 *      beta:     the final coef vector
 *      score:    the final score vector
 *      loglik:   the final log-likelihood
 *      n_iters:  number of iterations consumed
 *      error_status: 0 = no errors, -1 = input error, -2 = did not converge
 */
//' @title C++ function to fit interval censored models with elastic net
//'   penalty.
//'
//' @description \strong{This function is not meant to be called
//' on it's own!} Please use the \code{\link{iregnet}} function which includes
//' data validation and pre-processing.
//'
//' @param X Design matrix.
//' @param y Output matrix in a 2 column format. \code{NA}s denote censoring.
//' @param family String denoting the distribution to be fit.
//' @param lambda_path Vector containing the path of hyper-parameter lambda
//'   values.
//' @param debug Integer used for debugging during development. Unused otherwise.
//' @param out_status Vector containing censoring status of each observation.
//'   If not provided, it will be calculated using the \code{y} matrix.
//' @param intercept If \code{true}, intercept is to be fit. The first column
//'   of X must be \code{1}s if \code{intercept} is \code{true}.
//' @param alpha Hyper parameter for the elastic-net penalty.
//' @param scale_init The initial value of \code{scale} to be used for the fit.
//'   If not provided, a value is calculated depending on the distribution.
//' @param estimate_scale If \code{true}, \code{scale} is estimated. Else, the
//'     \code{scale} remains fixed at \code{scale_init}.
//' @param max_iter Maximum number of iterations to allow for \strong{each}
//'     model (\code{lambda} value) along the regularization path.
//' @param unreg_sol If \code{true}, the last model fit will be unregularized,
//'    i.e., the final \code{lambda} value will be \code{0}. Only used if
//'    \code{lambda_path} is not specified.
//' @param flag_standardize_x If \code{true}, the design matrix \code{X} will
//'     column-wise standardized.
//' @param threshold Convergence detected if absolute change in a coefficient
//'     is less than this value.
//' @param num_lambda Number of lambda values to use in the regularization
//'    path. Only used if \code{lambda_path} is not specified.
//' @param eps_lambda Ratio between the maximum and minimum values of
//'     \code{lambda}. Maximum value of \code{lambda} is calculated based on
//'     the distribution and the data. Only used if \code{lambda_path} is not
//'     specified.
//'
// [[Rcpp::export]]
Rcpp::List
fit_cpp(arma::mat& X, arma::mat& y,
        Rcpp::String family,   Rcpp::NumericVector lambda_path,
        int debug,             Rcpp::IntegerVector out_status,
        bool intercept,        double alpha,
        double scale_init,     bool estimate_scale,
        bool unreg_sol,        bool flag_standardize_x,
        double max_iter,       double threshold,
        int num_lambda,        double eps_lambda
        )
{
  /* Initialise some helper variables */
  IREG_DIST transformed_dist;  // Orig dist is the one that is initially given
                               // transformed_dist will be the dist of transformed output variables.
                               // Eg- orig_dist = "loglogistic", transformed_dist="logistic", with transform_y=log
                               // NOTE: This transformation is done before coming to this routine
  double scale;
  int error_status = 0;

  const ull n_obs  = X.n_rows;
  const ull n_vars = X.n_cols;  // n_vars is the number of variables corresponding to the coeffs of X
  transformed_dist = get_ireg_dist(family);

  /* Loggaussain = Gaussian with log(y), etc. */
  // get_transformed_dist(orig_dist, transformed_dist, &scale, &estimate_scale, y);

  /* Create output variables */
  if (lambda_path.size() > 0) {
    num_lambda = lambda_path.size();
  }
  Rcpp::NumericMatrix out_beta(n_vars, num_lambda + 1);       // will contain the entire series of solutions
  Rcpp::IntegerVector out_n_iters(num_lambda + 1);
  Rcpp::NumericVector out_lambda(num_lambda + 1);
  Rcpp::NumericVector out_scale(num_lambda + 1);
  Rcpp::NumericVector out_loglik(num_lambda + 1);

  /* use given values for the lambda path */
  if (lambda_path.size() > 0) {
    for(ull i = 0; i < num_lambda; ++i) {
      // Make sure that the given lambda_path is non-negative decreasing
      if (lambda_path[i] < 0 || (i > 0 && lambda_path[i] > lambda_path[i-1])) {
        Rcpp::stop("lambdas must be positive and decreasing.");
      }

      out_lambda[i] = lambda_path[i];
    }
  }

  double *beta = new double [n_vars];                           // Initially points to the first solution
  int *n_iters = INTEGER(out_n_iters);
  double *lambda_seq = REAL(out_lambda);

  double loglik = 0, scale_update = 0, log_scale;

  // TEMPORARY VARIABLES: not returned // TODO: Maybe alloc them together?
  rowvec eta_vec(n_obs, fill::zeros); // vector of linear predictors = X' beta
                                      // eta = 0 for the initial lambda_max, and in each iteration of coordinate descent,
                                      // eta is updated along with beta in place

  rowvec w_vec(n_obs, fill::zeros);   // diagonal of the hessian of LL wrt eta
                                      // these are the weights of the IRLS linear reg
  rowvec z_vec(n_obs, fill::zeros);   // z_i = eta_i - mu_i / w_i
  double *mu = new double [n_obs + 1];
  double *ym = new double [n_obs];    // the observation wise mean values for y's
  double mean_y;
  double *mean_x = new double [n_vars];
  double *std_x = new double [n_vars];
  IREG_CENSORING *status;
  int function_type = 1; // store the type of the function

  /* get censoring types of the observations */
  if (out_status.size() == 0) {
    status = new IREG_CENSORING [n_obs];
    get_censoring_types(y, status); // NANs denote censored observations in y, NOTE: y will be modified!
  }
  else {
    status = (IREG_CENSORING *) &out_status[0];
  }

  /* check out the whole obs by status[n_obs]
   * store the function type by censoring_types
   */
  for (int i = 0; i < n_obs; ++i) {
    if(status[i] == 3){
      function_type = 3;
      break;
    } else if(status[i] == 0 && function_type == 1) {
      function_type = 0;
    } else if(status[i] == 2 && function_type == 1){
      function_type = 2;
    }
  }

  /* X is columnwise variance normalized, mean is NOT set to 0
   * y is mean normalized.
   */
  /*if (flag_standardize_x) {
    standardize_x(X, mean_x, std_x, intercept);
  }*/
  mean_y = get_y_means(y, status, ym);
  standardize_y(y, ym, mean_y);

  /* SCALE RULES:
   * Whether or not it is estimated depends on estimate_scale. For exponential, this is forced to False and scale fixed to 1. // TODO
   *
   * If you provide no scale, a starting value will be calculated
   * If you provide a scale, it will be used as the initial value
   */
  if (scale_init == Rcpp::NA) {
    scale_init = get_init_var(ym, status, n_obs, transformed_dist);
    scale_init = 2 * sqrt(scale_init);    // use 2 * sqrt(var) as in survival
    // TODO: find corner cases where you need to take 2 * sqrt(var) as in survival
  }
  scale = scale_init; log_scale = log(scale);

  /* Sort the whole matrix by censoring type, get the sorted X and y
   */

  // Not support other distributions now
  if(transformed_dist != IREG_DIST_GAUSSIAN) {
    function_type = -1;
  }

  int *separator; // Store the index of separator in matrix X&y by censoring type.

  if(function_type != 1 && function_type != -1){

    mat sorted_X;
    mat sorted_y;

    separator = sort_X_and_y(sorted_X, sorted_y, status, X, y, function_type);

    X = sorted_X;
    y = sorted_y;
  }

  if (flag_standardize_x) {
    standardize_x(X, mean_x, std_x, intercept);
  }

  switch(function_type) {
    case 0:   target_compute_grad_response = compute_grad_response_gaussian_right;  break;
    case 1:   target_compute_grad_response = compute_grad_response_gaussian_none; break;
    case 2:   target_compute_grad_response = compute_grad_response_gaussian_left;   break;
    case 3:   target_compute_grad_response = compute_grad_response_gaussian_interval;   break;
    default:  target_compute_grad_response = compute_grad_response; break;
  }

  /* Initalize the pathwise solution
   * We will always start with lambda_max, at which beta = 0, eta = 0.
   * If some value of lambda is supplied, we will stop iterations at that value.
   */

  // set beta = 0 and eta = 0
  //beta = REAL(out_beta);
  for (ull i = 0; i < n_vars; ++i) {
    // sets only the first solution (first col of out_beta) to 0
    beta[i] = 0;
  }
  n_iters[0] = 0; out_scale[0] = scale;

  // Create Square X
  mat X_square = square(X);

  /******************************************************************************************/
  /* Iterate over grid of lambda values */
  bool flag_beta_converged = 0;
  rowvec sol_num_vec(n_vars, fill::zeros), sol_denom_vec(n_vars, fill::zeros);
  double beta_new;
  double old_scale;
  double lambda_max_unscaled;
  double eps_ratio = std::pow(eps_lambda, 1.0 / (num_lambda-1));
  // Separate Matrix y
  rowvec aram_y_l(n_obs);
  rowvec aram_y_r(n_obs);

  aram_y_l = (y.col(0)).t();

  if(y.n_cols==2){
    aram_y_r = (y.col(1)).t();
  }else{
    aram_y_r.resize(0);
  }

  // Optimize Temp Var
  rowvec w_division_nobs;
  vec flag_beta_converged_vec(n_vars, fill::zeros);// store converged result of beta[]
  rowvec temp_eta_vec(n_obs, fill::zeros);// store the temp result of (beta_new - beta_k) * X.col(k)
                                          // during COEFF LOOP
  mat eta_mat(n_obs, n_vars, fill::zeros);// save the result of every beta_k * X.col(k) in each column;
  vec temp_eta_vec_update(n_obs, fill::zeros);
  rowvec base_eta_vec;// store the result of (w/n_obs) * (z- eta) without (beta_new - beta_k) * X.col(k)
                      // during COEFF LOOP
  rowvec y_eta; // store the result of (y - eta)
  rowvec y_eta_square; // store the result of square(y - eta)

  /* below temp vars is used for distribution function
   * only support gaussian distribution & right censoring now
   */
  rowvec *compute_grad_response_temp_var = new rowvec [10];

  for (int m = 0; m < num_lambda + 1; ++m) {
    /* Compute the lambda path */
    if (lambda_path.size() == 0) {

      /* Do an initial fit with lambda set to BIG, will fit scale and intercept if applicable */
      if (m == 0) {
        lambda_max_unscaled = lambda_seq[0] = BIG;
      }

      /* Calculate lambda_max using intial scale fit */
      if (m == 1) {
        lambda_seq[m] = compute_lambda_max(X, &w_vec, &z_vec, &eta_vec, intercept, alpha, n_vars, n_obs, debug);

      /* Last solution should be unregularized if the flag is set */
      } else if (m == num_lambda && unreg_sol == true)
        lambda_seq[m] = 0;

      /* All other lambda calculated */
      else if (m > 1) {
        lambda_seq[m] = lambda_seq[m - 1] * eps_ratio;
      }
    }

    //This step is commented because the store way of beta is updated
    /* Initialize the solution at this lambda using previous lambda solution */
    // We need to explicitly do this because we need to store all the solutions separately
    /*if (m != 0) {                         // Initialise solutions using previous value
      for (ull i = 0; i < n_vars; ++i) {
        beta[i + n_vars] = beta[i];
      }
      beta = beta + n_vars;   // go to the next column
    }*/
    /* Before CYCLIC COORDINATE DESCENT, set the whole vector to zero*/
    flag_beta_converged_vec.zeros();
    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    n_iters[m] = 0;
    do {                                  // until Convergence of beta

      n_iters[m]++;

      flag_beta_converged = 1;            // = 1 if beta converges
      old_scale = scale;
      temp_eta_vec.zeros();

      y_eta = (aram_y_l - eta_vec) / scale;
      y_eta_square = square(y_eta);
      // IRLS: Reweighting step: calculate w and z again (beta & hence eta would have changed)  TODO: make dg ddg, local so that we can save computations?
      loglik = (*target_compute_grad_response)(&w_vec, &z_vec, &scale_update, &aram_y_l, &aram_y_r, eta_vec, scale,     // TODO:store a ptr to y?
                            status, n_obs, transformed_dist, NULL, debug==1 && m == 0, estimate_scale, y_eta, y_eta_square, separator, compute_grad_response_temp_var);

      // Calculate before iterate over beta elementwise loop
      w_division_nobs = w_vec / n_obs;
      sol_denom_vec = w_division_nobs * X_square;
      base_eta_vec = w_division_nobs % (z_vec - eta_vec);


      /* COEFF LOOP: iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_vars; ++k) {

        /*If current beta_k is already converged, continue to next beta iterated*/
        /*if(flag_beta_converged_vec(k) == 1)
          continue;
        */
        /* if eta_update_flag == 0, the eta is not changed, use the base_eta directly
         * if eta_update_flag == 1, the eta is changed, add the change part before calculate sol_num_vec(k)
         */
          sol_num_vec(k) = as_scalar(
                  (base_eta_vec + (w_division_nobs % ((eta_mat.col(k)).t() - temp_eta_vec)))
                  * X.col(k)
          );

        // Note: The signs given in the coxnet paper are incorrect, since the subdifferential should have a negative sign.
        sol_num_vec(k) *= -1; sol_denom_vec(k) *= -1;

        // if (debug == 1 && m == 0)
        //   std::cerr << n_iters[m] << " " << k << " " << "sols " << sol_num << " " << sol_denom << "\n";
        /* The intercept should not be regularized, and hence is calculated directly */
        if (intercept && k == 0) {
          beta_new = sol_num_vec(k) / sol_denom_vec(k);

        } else {
          beta_new = soft_threshold(sol_num_vec(k), lambda_seq[m] * alpha) /
                     (sol_denom_vec(k) + lambda_seq[m] * (1 - alpha));
        }

        // if any beta_k has not converged, we will come back for another cycle.
	double abs_change = fabs(beta_new - beta[k]);
        if (abs_change > threshold) {
	  if(debug==1 && max_iter==n_iters[m])printf("iter=%d lambda=%d beta_%lld not converged, abs_change=%f > %f=threshold\n", n_iters[m], m, k, abs_change, threshold);
          flag_beta_converged = 0;
          temp_eta_vec_update = (beta_new - beta[k]) * X.col(k);
          temp_eta_vec += temp_eta_vec_update.t();// this will contain the new beta_k
          eta_mat.col(k) += temp_eta_vec_update;// save the result of beta_k * X.col(k)
          beta[k] = beta_new;
        } else {
          /*The current beta_k has converged during this CYCLIC COORDINATE DESCENT*/
          flag_beta_converged_vec(k) = 1;
        }

        // if (debug==1 && m == 1)
        //   std::cerr << n_iters[m] << " " << k << " " << " BETA " << beta[k] << "\n";

        //eta[] is already updated in the below step
        /*for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
          // if (debug==1 && m==0) {
          //   std::cerr << n_iters[m] << " " << i << " " << "ETA" <<  eta[i] << "\n";
          // }
        }*/
      }   // end for: beta_k solution
      eta_vec += temp_eta_vec;// Add the total change to eta for next computation

      if (estimate_scale) {
        log_scale += scale_update; scale = exp(log_scale);

        // if (fabs(scale - old_scale) > threshold) {    // TODO: Maybe should be different for sigma?
	double abs_change = fabs(scale - old_scale);
        if (abs_change > threshold) {    // TODO: Maybe should be different for sigma?
	  if(debug==1 && max_iter==n_iters[m])printf("iter=%d lambda=%d scale not converged, abs_change=%f > %f=threshold\n", n_iters[m], m, abs_change, threshold);
          flag_beta_converged = 0;
        }
      }
    } while ((flag_beta_converged != 1) && (n_iters[m] < max_iter));

    out_loglik[m] = loglik;
    out_scale[m] = scale;

    /* Check for errors */
    if (n_iters[m] == max_iter)
      error_status = -1;
    if (std::isinf(out_loglik[m]))
      error_status = -3;
    if (std::isnan(out_loglik[m])) {  // Fatal error: If NaNs are produced something is wrong.
      error_status = 1;
      Rcpp::stop("NANs produced");
    }
    //Store new beta[]
    std::copy(beta, beta + n_vars, REAL(out_beta) + (m) * n_vars);
  } // end for: lambda

  /* Scale the coefs back to the original scale */
  for (ull m = 0; m < num_lambda + 1; ++m) {
    //if (transformed_dist == IREG_DIST_LOGISTIC)
      out_beta(0, m) += mean_y;     // intercept will contain the contribution of mean_y

    if (flag_standardize_x) {
      double intercept_diff = 0;
      for (ull i = int(intercept); i < n_vars; ++i) {  // no standardization for intercept
        out_beta(i, m) = out_beta(i, m) / std_x[i];
        intercept_diff += out_beta(i, m) * mean_x[i];
      }
      out_beta(0, m) -= intercept_diff;
    }
  }

  /* Free the temporary variables */
  delete [] beta;
  delete [] mu;
  delete [] ym;
  delete [] compute_grad_response_temp_var;
  if (out_status == Rcpp::NA) // we would've allocated a new vector in this case
    delete [] status;

  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta(Rcpp::_, Rcpp::Range(1,num_lambda)),
                            Rcpp::Named("lambda")       = out_lambda[Rcpp::Range(1,num_lambda)],
                            Rcpp::Named("num_lambda")   = num_lambda,
                            Rcpp::Named("n_iters")      = out_n_iters[Rcpp::Range(1,num_lambda)],
                            Rcpp::Named("loglik")       = out_loglik[Rcpp::Range(1,num_lambda)],
                            Rcpp::Named("scale")        = out_scale[Rcpp::Range(1,num_lambda)],
                            Rcpp::Named("estimate_scale") = estimate_scale,
                            Rcpp::Named("scale_init")   = scale_init,
                            Rcpp::Named("error_status") = error_status
                            );
}

static inline double
get_y_means (mat &y, IREG_CENSORING *status, double *ym)
{
  if (!ym || !status)
    return 0;

  double mean_y = 0;
  for (ull i = 0; i < y.n_rows; ++i) {
      switch (status[i]) {
        case IREG_CENSOR_LEFT:
        case IREG_CENSOR_RIGHT:
        case IREG_CENSOR_NONE:
          ym[i] = y(i, 0);
          break;

        case IREG_CENSOR_INTERVAL:
          ym[i] = (y(i, 1) + y(i, 0)) / 2;
          break;

        default:
          ym[i] = 0;
      }
      mean_y += ym[i];
  } // end for
  mean_y /= y.n_rows;
  return mean_y;
}

void
get_censoring_types (mat &y, IREG_CENSORING *status)
{
  double y_l, y_r;
  for (ull i = 0; i < y.n_rows; ++i) {
    if (std::isinf(fabs(y(i, 0)))) y(i, 0) = NAN;
    if (std::isinf(fabs(y(i, 1)))) y(i, 1) = NAN;
    y_l = y(i, 0); y_r = y(i ,1);

    if (y_l == Rcpp::NA) {
      if (y_r == Rcpp::NA)
        Rcpp::stop("Invalid interval: both limits NA");
      else {
        status[i] = IREG_CENSOR_LEFT;       // left censoring
        y(i, 0) = y_r; // NOTE: We are putting the value in the left col, survival style!
      }
      continue;
    }

    if (y_r == Rcpp::NA)                // right censoring
      status[i] = IREG_CENSOR_RIGHT;
    else if (y_l == y_r)
      status[i] = IREG_CENSOR_NONE;         // no censoring
    else if (y_l > y_r)
      Rcpp::stop("Invalid interval: start > stop");
    else
      status[i] = IREG_CENSOR_INTERVAL;     // interval censoring

  } // end for
}


/*
 * Takes as input Rcpp string family name, and returns corresponding enum val
 */
IREG_DIST
get_ireg_dist (Rcpp::String dist_str)
{
  if (strcmp("gaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_GAUSSIAN;
  if (strcmp("logistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOGISTIC;
  if (strcmp("extreme_value", dist_str.get_cstring()) == 0)
    return IREG_DIST_EXTREME_VALUE;
  if (strcmp("exponential", dist_str.get_cstring()) == 0)
    return IREG_DIST_EXPONENTIAL;
  if (strcmp("weibull", dist_str.get_cstring()) == 0)
    return IREG_DIST_WEIBULL;
  if (strcmp("loggaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOG_GAUSSIAN;
  if (strcmp("loglogistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOG_LOGISTIC;
  return IREG_DIST_UNKNOWN;
}

static inline double
soft_threshold (double x, double lambda)
{
  double temp = fabs(x) - lambda;
  return (temp > 0)? ((x > 0)? 1: -1) * temp: 0;
}

static void
standardize_x (mat &X,
               double *mean_x, double *std_x,
               bool intercept)
{
  for (ull i = int(intercept); i < X.n_cols; ++i) {  // don't standardize intercept col.
    mean_x[i] = std_x[i] = 0;
    for (ull j = 0; j < X.n_rows; ++j) {
      mean_x[i] += X(j, i);
    }
    mean_x[i] /= X.n_rows;

    for (ull j = 0; j < X.n_rows; ++j) {
      X(j, i) = X(j, i) - mean_x[i];  // center X
      std_x[i] += X(j, i) * X(j, i);
    }

    // not using (N-1) in denominator to agree with glmnet
    std_x[i] = sqrt(std_x[i] / X.n_rows);

    for (ull j = 0; j < X.n_rows; ++j) {
      X(j, i) = X(j, i) / std_x[i];
    }
  }
}

static void
standardize_y (mat &y, double *ym, double &mean_y)
{
  for (ull i = 0; i < y.n_rows * y.n_cols; ++i) {
    y[i] -= mean_y;
  }
  for (ull i = 0; i < y.n_rows; ++i)
    ym[i] -= mean_y;
}

static double
get_init_var (double *ym, IREG_CENSORING *status, ull n, IREG_DIST dist)
{
  double mean = 0, var = 0;

  for (int i = 0; i < n; ++i) {
    mean += ym[i];
  }
  mean = mean / n;

  for (int i = 0; i < n; ++i) {
    var += pow((ym[i] - mean), 2);
  }
  var = var / n;

  switch (dist) {
    case IREG_DIST_EXTREME_VALUE :
      mean = mean + 0.572;
      var = var / 1.64;
    break;

    case IREG_DIST_LOGISTIC:
      var = var / 3.2;
    break;
  }

  return var;
}

static inline double
compute_lambda_max(mat X, rowvec *w, rowvec *z, rowvec *eta,
                   bool intercept, double &alpha, ull n_vars, ull n_obs,
                   bool debug=0)
{
  double lambda_max = -1; // calculated values are always non-negative
  for (ull j = int(intercept); j < n_vars; ++j) {   // dont include intercept col in lambda calc.
    double temp = 0;

    for (ull i = 0; i < n_obs; ++i) {
      // NOTE: `eta` contains only the contribution of the intercept, since all
      // other `beta` values are 0.
      temp += ((*w)(i) * X(i, j) * ((*z)(i) - (*eta)(i)));
    }
    temp = fabs(temp);
    // if (debug) {
    //   std::cerr << "LAMBDA " << temp / n_obs << "\n";
    // }
    lambda_max = (lambda_max > temp)? lambda_max: temp;
  }
  lambda_max /= (n_obs * max(alpha, 1e-3));  // prevent divide by zero

  return lambda_max;
}

static int *
sort_X_and_y(mat &sorted_X, mat &sorted_y, IREG_CENSORING *status, mat X, mat y, int function_type)
{
  int *separator = new int [4]; // Store the index of separator in matrix X&y by censoring type.

  if(function_type != 3){

    /*
     * Temp var for sort the matrix X&y
     */
    mat none_censoring_X_mat = mat(X.n_rows, X.n_cols);
    mat left_or_right_censoring_X_mat = mat(X.n_rows, X.n_cols);
    mat none_censoring_y_mat = mat(y.n_rows, y.n_cols);
    mat left_or_right_censoring_y_mat = mat(y.n_rows, y.n_cols);

    ull none_censoring_number = 0;
    ull left_or_right_censoring_number = 0;

    /*
     * Sort by censoring type,
     */
    for (ull i = 0; i < X.n_rows; ++i) {

      if(status[i] == 1){

//        none_censoring_X_mat.insert_rows(none_censoring_number, X.row(i));
//        none_censoring_y_mat.insert_rows(none_censoring_number, y.row(i));
        none_censoring_X_mat.row(none_censoring_number) = X.row(i);
        none_censoring_y_mat.row(none_censoring_number) = y.row(i);

        none_censoring_number++;
      } else {

//        left_or_right_censoring_X_mat.insert_rows(left_or_right_censoring_number, X.row(i));
//        left_or_right_censoring_y_mat.insert_rows(left_or_right_censoring_number, y.row(i));

        left_or_right_censoring_X_mat.row(left_or_right_censoring_number) = X.row(i);
        left_or_right_censoring_y_mat.row(left_or_right_censoring_number) = y.row(i);

        left_or_right_censoring_number++;
      }
    }

    none_censoring_X_mat.resize(none_censoring_number, X.n_cols);
    none_censoring_y_mat.resize(none_censoring_number, y.n_cols);

    left_or_right_censoring_X_mat.resize(left_or_right_censoring_number, X.n_cols);
    left_or_right_censoring_y_mat.resize(left_or_right_censoring_number, y.n_cols);

    none_censoring_X_mat.insert_rows(none_censoring_number, left_or_right_censoring_X_mat);
    none_censoring_y_mat.insert_rows(none_censoring_number, left_or_right_censoring_y_mat);

    sorted_X = none_censoring_X_mat;
    sorted_y = none_censoring_y_mat;

    separator[0] = none_censoring_number;

  } else {

    /*
     * Temp var for sort the matrix X&y
     */
    mat none_censoring_X_mat = mat(X.n_rows, X.n_cols);
    mat left_censoring_X_mat = mat(X.n_rows, X.n_cols);
    mat right_censoring_X_mat = mat(X.n_rows, X.n_cols);
    mat interval_censoring_X_mat = mat(X.n_rows, X.n_cols);

    mat none_censoring_y_mat = mat(y.n_rows, y.n_cols);
    mat left_censoring_y_mat = mat(y.n_rows, y.n_cols);
    mat right_censoring_y_mat = mat(y.n_rows, y.n_cols);
    mat interval_censoring_y_mat = mat(y.n_rows, y.n_cols);

    ull none_censoring_number = 0;
    ull left_censoring_number = 0;
    ull right_censoring_number = 0;
    ull interval_censoring_number = 0;

    /*
     * Sort by censoring type,
     */
    for (ull i = 0; i < X.n_rows; ++i) {

      if(status[i] == 1){

        none_censoring_X_mat.row(none_censoring_number) = X.row(i);
        none_censoring_y_mat.row(none_censoring_number) = y.row(i);

        none_censoring_number++;
      } else if(status[i] == 0){

        right_censoring_X_mat.row(right_censoring_number) = X.row(i);
        right_censoring_y_mat.row(right_censoring_number) = y.row(i);

        right_censoring_number++;
      } else if(status[i] == 2){

        left_censoring_X_mat.row(left_censoring_number) = X.row(i);
        left_censoring_y_mat.row(left_censoring_number) = y.row(i);

        left_censoring_number++;
      } else {

        interval_censoring_X_mat.row(interval_censoring_number) = X.row(i);
        interval_censoring_y_mat.row(interval_censoring_number) = y.row(i);

        interval_censoring_number++;
      }
    }

    none_censoring_X_mat.resize(none_censoring_number, X.n_cols);
    none_censoring_y_mat.resize(none_censoring_number, y.n_cols);

    left_censoring_X_mat.resize(left_censoring_number, X.n_cols);
    left_censoring_y_mat.resize(left_censoring_number, y.n_cols);

    right_censoring_X_mat.resize(right_censoring_number, X.n_cols);
    right_censoring_y_mat.resize(right_censoring_number, y.n_cols);

    interval_censoring_X_mat.resize(interval_censoring_number, X.n_cols);
    interval_censoring_y_mat.resize(interval_censoring_number, y.n_cols);

    none_censoring_X_mat.insert_rows(none_censoring_number, left_censoring_X_mat);
    none_censoring_X_mat.insert_rows(none_censoring_number + left_censoring_number, right_censoring_X_mat);
    none_censoring_X_mat.insert_rows(none_censoring_number + left_censoring_number + right_censoring_number,
                                     interval_censoring_X_mat);

    none_censoring_y_mat.insert_rows(none_censoring_number, left_censoring_y_mat);
    none_censoring_y_mat.insert_rows(none_censoring_number + left_censoring_number, right_censoring_y_mat);
    none_censoring_y_mat.insert_rows(none_censoring_number + left_censoring_number + right_censoring_number,
                                     interval_censoring_y_mat);

    sorted_X = none_censoring_X_mat;
    sorted_y = none_censoring_y_mat;

    separator[0] = none_censoring_number;
    separator[1] = left_censoring_number;
    separator[2] = right_censoring_number;
    separator[3] = interval_censoring_number;
  }

  return separator;
}