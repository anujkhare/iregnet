/* for easy compilation in emacs -*- compile-command: "R CMD INSTALL .." -*- */

/*
 * iregnet_fit.cpp
 * Author: Anuj Khare <khareanuj18@gmail.com>
 * We use 2 spaces per tab, and expand the tabs
 */
#include "iregnet.h"
#include <R_ext/Utils.h>

#define BIG 1e35

static inline double
get_y_means (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *ym);
static inline double soft_threshold (double x, double lambda);
static double get_init_var (double *ym, IREG_CENSORING *status, ull n, IREG_DIST dist);
static void standardize_x (Rcpp::NumericMatrix &X,
                           double *mean_x, double *std_x,
                           bool intercept);
static void
standardize_y (Rcpp::NumericMatrix &y, double *ym, double &mean_y);
static inline double compute_lambda_max(Rcpp::NumericMatrix X, double *w, double *z,
                                        double *eta, bool intercept, double &alpha,
                                        ull n_vars, ull n_obs, bool debug);

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
fit_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y,
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
  bool flag_lambda_given = (lambda_path.size() > 0);
  int error_status = 0;

  const ull n_obs  = X.nrow();
  const ull n_vars = X.ncol();  // n_vars is the number of variables corresponding to the coeffs of X
  transformed_dist = get_ireg_dist(family);

  /* Loggaussain = Gaussian with log(y), etc. */
  // get_transformed_dist(orig_dist, transformed_dist, &scale, &estimate_scale, y);

  /* Create output variables */
  if (flag_lambda_given) {
    num_lambda = lambda_path.size();
  }
  Rcpp::NumericMatrix out_beta(n_vars, num_lambda + 1);       // will contain the entire series of solutions
  Rcpp::IntegerVector out_n_iters(num_lambda + 1);
  Rcpp::NumericVector out_lambda(num_lambda + 1);
  Rcpp::NumericVector out_scale(num_lambda + 1);
  Rcpp::NumericVector out_loglik(num_lambda + 1);

  /* use given values for the lambda path */
  if (flag_lambda_given) {
    for(ull i = 0; i < num_lambda; ++i) {
      // Make sure that the given lambda_path is non-negative decreasing
      if (lambda_path[i] < 0 || (i > 0 && lambda_path[i] > lambda_path[i-1])) {
        Rcpp::stop("lambdas must be positive and decreasing.");
      }

      out_lambda[i] = lambda_path[i];
    }
  }

  double *beta;                           // Initially points to the first solution
  int *n_iters = INTEGER(out_n_iters);
  double *lambda_seq = REAL(out_lambda);

  double loglik = 0, scale_update = 0, log_scale;

  // TEMPORARY VARIABLES: not returned // TODO: Maybe alloc them together?
  double *eta = new double [n_obs];   // vector of linear predictors = X' beta
                                      // eta = 0 for the initial lambda_max, and in each iteration of coordinate descent,
                                      // eta is updated along with beta in place

  double *w  = new double [n_obs];    // diagonal of the hessian of LL wrt eta
                                      // these are the weights of the IRLS linear reg
  double *z = new double [n_obs];     // z_i = eta_i - mu_i / w_i
  double *mu = new double [n_obs + 1];
  double *ym = new double [n_obs];    // the observation wise mean values for y's
  double mean_y;
  double *mean_x = new double [n_vars];
  double *std_x = new double [n_vars];
  IREG_CENSORING *status;


  /* get censoring types of the observations */
  if (out_status.size() == 0) {
    status = new IREG_CENSORING [n_obs];
    get_censoring_types(y, status); // NANs denote censored observations in y, NOTE: y will be modified!
  }
  else {
    status = (IREG_CENSORING *) &out_status[0];
  }

  /* X is columnwise variance normalized, mean is NOT set to 0
   * y is mean normalized.
   */
  if (flag_standardize_x) {
    standardize_x(X, mean_x, std_x, intercept);
  }
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

  /* Initalize the pathwise solution
   * We will always start with lambda_max, at which beta = 0, eta = 0.
   * If some value of lambda is supplied, we will stop iterations at that value.
   */

  // set beta = 0 and eta = 0
  beta = REAL(out_beta);
  for (ull i = 0; i < n_vars; ++i) {
    // sets only the first solution (first col of out_beta) to 0
    beta[i] = 0;
  }
  n_iters[0] = 0; out_scale[0] = scale;

  for (ull i = 0; i < n_obs; ++i) {
    eta[i] = w[i] = z[i] = 0;
  }

  /******************************************************************************************/
  /* Iterate over grid of lambda values */
  bool flag_beta_converged = 0;
  double sol_num, sol_denom;
  double beta_new;
  double old_scale;
  double lambda_max_unscaled;
  double eps_ratio = std::pow(eps_lambda, 1.0 / (num_lambda-1));
  // There is an extra lambda for initial_fit if we calculate them.
  int end_ind = num_lambda + !flag_lambda_given - 1;

  for (int m = 0; m <= end_ind; ++m) {
    /* Compute the lambda path */
    if (!flag_lambda_given) {

      /* Do an initial fit with lambda set to BIG, will fit scale and intercept if applicable */
      if (m == 0) {
        lambda_max_unscaled = lambda_seq[0] = BIG;
      }

      /* Calculate lambda_max using intial scale fit */
      if (m == 1) {
        lambda_seq[m] = compute_lambda_max(X, w, z, eta, intercept, alpha, n_vars, n_obs, debug);

      /* Last solution should be unregularized if the flag is set */
      } else if (m == end_ind && unreg_sol == true)
        lambda_seq[m] = 0;

      /* All other lambda calculated */
      else if (m > 1) {
        lambda_seq[m] = lambda_seq[m - 1] * eps_ratio;
      }
    }

    /* Initialize the solution at this lambda using previous lambda solution */
    // We need to explicitly do this because we need to store all the solutions separately
    if (m != 0) {                         // Initialise solutions using previous value
      for (ull i = 0; i < n_vars; ++i) {
        beta[i + n_vars] = beta[i];
      }
      beta = beta + n_vars;   // go to the next column
    }

    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    n_iters[m] = 0;
    do {                                  // until Convergence of beta
      R_CheckUserInterrupt(); // Ctrl-C to interrupt.
      
      n_iters[m]++;

      flag_beta_converged = 1;            // = 1 if beta converges
      old_scale = scale;

      // IRLS: Reweighting step: calculate w and z again (beta & hence eta would have changed)  TODO: make dg ddg, local so that we can save computations?
      // ***Vectorize
      loglik = compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                            status, n_obs, transformed_dist, NULL, debug==1 && m == 0);
      /* iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_vars; ++k) {
        sol_num = sol_denom = 0;
        // ***Equation 18
        // ***Vectorize
        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
          sol_num += (w[i] * X(i, k) * (z[i] - eta[i])) / n_obs;
          sol_denom += (w[i] * X(i, k) * X(i, k)) / n_obs;
        }

        // Note: The signs given in the coxnet paper are incorrect, since the subdifferential should have a negative sign.
        sol_num *= -1; sol_denom *= -1;

        // if (debug == 1 && m == 0)
        //   std::cerr << n_iters[m] << " " << k << " " << "sols " << sol_num << " " << sol_denom << "\n";
        /* The intercept should not be regularized, and hence is calculated directly */
        if (intercept && k == 0) {
          beta_new = sol_num / sol_denom;

        } else { // ***Eqn 19
          beta_new = soft_threshold(sol_num, lambda_seq[m] * alpha) /
                     (sol_denom + lambda_seq[m] * (1 - alpha));
        }

        // if any beta_k has not converged, we will come back for another cycle.
	double abs_change = fabs(beta_new - beta[k]);
        if (abs_change > threshold) {
	  if(debug==1 && max_iter==n_iters[m])printf("iter=%d lambda=%d beta_%lld not converged, abs_change=%f > %f=threshold\n", n_iters[m], m, k, abs_change, threshold);
          flag_beta_converged = 0;
          beta[k] = beta_new;
        }

        // if (debug==1 && m == 1)
        //   std::cerr << n_iters[m] << " " << k << " " << " BETA " << beta[k] << "\n";
        // ***Vectorize
        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
          // if (debug==1 && m==0) {
          //   std::cerr << n_iters[m] << " " << i << " " << "ETA" <<  eta[i] << "\n";
          // }
        }

      }   // end for: beta_k solution

      if (estimate_scale) {
        log_scale += scale_update; scale = exp(log_scale);

        // if (fabs(scale - old_scale) > threshold) {    // TODO: Maybe should be different for sigma?
	        double abs_change = fabs(scale - old_scale);
          if (abs_change > threshold) {    // TODO: Maybe should be different for sigma?
	          if(debug==1 && max_iter==n_iters[m])printf("iter=%d lambda=%d scale not converged, abs_change=%f > %f=threshold\n", n_iters[m], m, abs_change, threshold);
              flag_beta_converged = 0;
        }
      }
      // flag_beta_converged = 1 then converged
    } while ((flag_beta_converged != 1) && (n_iters[m] < max_iter));

    out_loglik[m] = loglik;
    out_scale[m] = scale;

    /* Check for errors */
    if (n_iters[m] == max_iter)
      error_status = -1;
    if (std::isnan(out_loglik[m])) {  // Fatal error: If NaNs are produced something is wrong.
      error_status = -3;
      break;
    }
  } // end for: lambda

  if(error_status == -1)
    Rcpp::warning("Ran out of iterations and failed to converge.");
  if(error_status == -3)
   Rcpp::warning("Failed to converge. Try again after adding more data.");

  /* Scale the coefs back to the original scale */
  for (ull m = 0; m <= end_ind; ++m) {
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
  delete [] eta;
  delete [] mu;
  delete [] w;
  delete [] z;
  delete [] ym;
  if (out_status == Rcpp::NA) // we would've allocated a new vector in this case
    delete [] status;

  int start_ind = !flag_lambda_given;
  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta(Rcpp::_, Rcpp::Range(start_ind ,end_ind)),
                            Rcpp::Named("lambda")       = out_lambda[Rcpp::Range(start_ind ,end_ind)],
                            Rcpp::Named("num_lambda")   = num_lambda,
                            Rcpp::Named("n_iters")      = out_n_iters[Rcpp::Range(start_ind ,end_ind)],
                            Rcpp::Named("loglik")       = out_loglik[Rcpp::Range(start_ind ,end_ind)],
                            Rcpp::Named("scale")        = out_scale[Rcpp::Range(start_ind ,end_ind)],
                            Rcpp::Named("estimate_scale") = estimate_scale,
                            Rcpp::Named("scale_init")   = scale_init,
                            Rcpp::Named("error_status") = error_status
                            );
}

static inline double
get_y_means (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *ym)
{
  if (!ym || !status)
    return 0;

  double mean_y = 0;
  for (ull i = 0; i < y.nrow(); ++i) {
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
  mean_y /= y.nrow();
  return mean_y;
}

void
get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *status)
{
  double y_l, y_r;
  for (ull i = 0; i < y.nrow(); ++i) {
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
standardize_x (Rcpp::NumericMatrix &X,
               double *mean_x, double *std_x,
               bool intercept)
{
  for (ull i = int(intercept); i < X.ncol(); ++i) {  // don't standardize intercept col.
    mean_x[i] = std_x[i] = 0;
    for (ull j = 0; j < X.nrow(); ++j) {
      mean_x[i] += X(j, i);
    }
    mean_x[i] /= X.nrow();

    for (ull j = 0; j < X.nrow(); ++j) {
      X(j, i) = X(j, i) - mean_x[i];  // center X
      std_x[i] += X(j, i) * X(j, i);
    }

    // not using (N-1) in denominator to agree with glmnet
    std_x[i] = sqrt(std_x[i] / X.nrow());

    for (ull j = 0; j < X.nrow(); ++j) {
      X(j, i) = X(j, i) / std_x[i];
    }
  }
}

static void
standardize_y (Rcpp::NumericMatrix &y, double *ym, double &mean_y)
{
  for (ull i = 0; i < y.nrow() * y.ncol(); ++i) {
    y[i] -= mean_y;
  }
  for (ull i = 0; i < y.nrow(); ++i)
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
compute_lambda_max(Rcpp::NumericMatrix X, double *w, double *z, double *eta,
                   bool intercept, double &alpha, ull n_vars, ull n_obs,
                   bool debug=0)
{
  double lambda_max = -1; // calculated values are always non-negative
  for (ull j = int(intercept); j < n_vars; ++j) {   // dont include intercept col in lambda calc.
    double temp = 0;

    for (ull i = 0; i < n_obs; ++i) {
      // NOTE: `eta` contains only the contribution of the intercept, since all
      // other `beta` values are 0.
      temp += (w[i] * X(i, j) * (z[i] - eta[i]));
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
