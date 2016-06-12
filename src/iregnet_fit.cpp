#include "iregnet.h"

#define BIG 1e35

static inline void get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *censoring_type);
static inline double soft_threshold(double x, double lambda);
static void standardize_x_y(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                            double *mean_x, double *std_x, double &mean_y,
                            double &std_y, bool intercept);

double identity (double y)
{
  return y;
}

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
                   bool intercept,
                   double scale,     // bool estimate_scale = false,    // TODO: SCALE??
                bool standardize = false,   
		double max_iter = 1000,  double threshold = 1e-4,
                   int num_lambda = 100,  double eps_lambda = 0.0001,   // TODO: depends on nvars vs nobs
                   int flag_debug = 0)
{
  bool estimate_scale = 0;
  if (scale == Rcpp::NA) {
    std::cout << "Estimating scale\n";
    estimate_scale = 1;
    // TODO: NEED INITIAL VALUE OF SCALE NOW!
    scale = 10;
  }

  /* Uselesss print stuff for now */
  if (flag_debug == IREG_DEBUG_INPUT) {       // 1
    // Rcpp has implemented << for NumericVectors and matrices too!
    std::cout << "X:\n" << X << std::endl;
    std::cout << "y:\n" << y << std::endl;
    std::cout << "family:" << family.get_cstring() << std::endl;
    std::cout << "alpha:" << alpha << std::endl;
  }


  /* Initialise some helper variables */
  ull n_obs, n_vars, n_cols_x, n_cols_y, n_params;
  IREG_DIST orig_dist, transformed_dist;  // Orig dist is the one that is initially given
                                          // transformed_dist will be the dist of transformed output variables. Eg- orig_dist = "loglogistic", transformed_dist="logistic", with transform_y=log
  double (*transform_y) (double y);

  n_obs  = X.nrow();
  n_cols_x = X.ncol();
  n_cols_y = y.ncol();
  orig_dist = get_ireg_dist(family);

  n_vars = X.ncol();  // n_vars is the number of variables corresponding to the coeffs of X
  n_params = n_vars + int(estimate_scale); // n_params is the number of parameters
                                           // to optimize (includes scale as well)

  if (flag_debug == IREG_DEBUG_N) {         // 2
    std::cout << "n_vars: " << n_vars << ", n_params: " << n_params << "\n";
    std::cout << "n_obs: " << n_obs << "\n";
  }


  /* TODO: Validate all the arguments again */
  if ((alpha > 1 || alpha < 0) || (IREG_DIST_UNKNOWN == orig_dist) ||
      (y.nrow() != n_obs)) {
    return Rcpp::List::create(Rcpp::Named("error_status") = -1);
  }


  /*
   * TODO: a function to apply distribution specific transformations and parameters
   * should preceed this fit function.
   */
  /* Apply the required transformation to the time scale (t = log(y))
   * the transformation depends on the distribution used.
   */
  switch(orig_dist) {
    case IREG_DIST_EXTREME_VALUE:
    case IREG_DIST_GAUSSIAN:
    case IREG_DIST_LOGISTIC:
      transform_y = identity;
      transformed_dist = orig_dist;
      break;

    case IREG_DIST_LOG_GAUSSIAN:
      transform_y = log;
      transformed_dist = IREG_DIST_GAUSSIAN;
      break;

    case IREG_DIST_EXPONENTIAL:
      transform_y = log;
      transformed_dist = IREG_DIST_EXTREME_VALUE;
      scale = 1;
      estimate_scale = 0;
      break;
  }

  // TODO: unnecessary computations
  for (ull i = 0; i < 2 * n_obs; ++i) {     // Rcpp::NumericMatrix is unrolled col major
    y[i] = transform_y(y[i]);
  }

  /*
   * Standardize functions:
   * Normalize X and y
   * TODO: Should you do that? How does it work for censored data?
   */
  double *mean_x = new double [n_vars], mean_y;
  double *std_x = new double [n_vars], std_y;

  /* NOTE:
   * we have scaled y and x, so you need to scale the obtained lambda values,
   * and coef (beta) values back to the original scale before returning them.
   */
  if (standardize) {
    standardize_x_y(X, y, mean_x, std_x, mean_y, std_y, intercept);      // FIXME: so that values are always estimated as for intercepts
    // standardize_x_y(X, y, mean_x, std_x, mean_y, std_y, true);
    std::cout << "y\n" << y << "\nx\n" << X;
    std::cout << "mean_y " << mean_y << "std_y " << std_y << "mean_x:\n";
    for (ull i=0; i<n_vars; ++i) {
      std::cout << i << " " << mean_x[i] << " " << std_x[i] << "\n";
    }
  }

  /* Create output variables */
  // Rcpp::NumericVector out_beta(n_params);
  Rcpp::NumericMatrix out_beta(n_params, num_lambda + 1);       // will contain the entire series of solutions
  Rcpp::IntegerVector out_n_iters(num_lambda + 1);
  Rcpp::NumericVector out_lambda(num_lambda + 1);
  Rcpp::NumericVector out_scale(num_lambda + 1);
  //Rcpp::NumericVector out_intercept(num_lambda, 0.0);         // may contain non-zero values only if intercept==T
==== BASE ====

  double *beta  = REAL(out_beta);                           // Initially points to the first solution
  int *n_iters = INTEGER(out_n_iters);
  double *lambda_seq = REAL(out_lambda);

  double loglik = 0, scale_update = 0, log_scale = log(scale);


  // TEMPORARY VARIABLES: not returned // TODO: Maybe alloc them together?
  double *eta = new double [n_obs];   // vector of linear predictors = X' beta
                                      // eta = 0 for the initial lambda_max, and in each iteration of coordinate descent,
                                      // eta is updated along with beta in place

  double *w  = new double [n_obs];    // diagonal of the hessian of LL wrt eta
                                      // these are the weights of the IRLS linear reg
  double *z = new double [n_obs];     // z_i = eta_i - mu_i / w_i

  /* get censoring types of the observations */ /* TODO: Incorporate survival style censoring */
  IREG_CENSORING censoring_type[n_obs];
  get_censoring_types(y, censoring_type); // NANs denote censored observations in y

  if (flag_debug == IREG_DEBUG_CENSORING) {         // 3
    Rcpp::NumericVector rvec (censoring_type, censoring_type + n_obs);
    std::cout << rvec << "\n";
    return Rcpp::List::create(Rcpp::Named("error_status")   = 0,
                              Rcpp::Named("censoring_types") = rvec);
  }


  /* Initalize the pathwise solution
   * We will always start with lambda_max, at which beta = 0, eta = 0.
   * If some value of lambda is supplied, we will stop iterations at that value.
   */

  // set beta = 0 and eta = 0
  for (ull i = 0; i < n_params; ++i) {
    // sets only the first solution (first col of out_beta) to 0
    beta[i] = 0;
  }
  n_iters[0] = 0; out_scale[0] = scale;

  for (ull i = 0; i < n_obs; ++i) {
    eta[i] = w[i] = z[i] = 0;
  }

  /////////////////////////////////////////
  // Calculate w and z right here!
  compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,
                        censoring_type, n_obs, transformed_dist, NULL);
  // std::cout <<" w z \n";
  // for (ull i=0; i<n_obs; ++i) {
  //   std::cout << w[i] << " " << z[i] << "\n";
  // }
  // std::cout << "\n";
  // std::cout << "\n";

  // Calculate lambda_max
    // TODO: try to optimize by reversing j and i loops and using an extra array
  // First, start with lambda_max = BIG (really really big), so that eta and beta are surely 0
  lambda_seq[0] = -BIG;
  for (ull j = int(intercept); j < n_vars; ++j) {   // dont include intercept col in lambda calc.
    double temp = 0;
    double tt = 0;

    for (ull i = 0; i < n_obs; ++i) {
      temp += (w[i] * X(i, j) * z[i]);

      tt += (X(i,j) * y[i]);
    }
    //temp = fabs(temp) / (n_obs * alpha);
    temp = fabs(temp) / (alpha);              // FIXME: It seems that is how they have calc lambda in GLMNET
    // std::cout << temp << " " << n_obs << " " << tt << "\n";

    lambda_seq[0] = (lambda_seq[0] > temp)? lambda_seq[0]: temp;
  }

  // TODO: you should only set lambda_max = Inf, and let it calc beta = eta = 0 itself.

  /* Iterate over grid of lambda values */
  double eps_ratio = std::pow(eps_lambda, 1.0 / (num_lambda-1));
  bool flag_beta_converged = 0;
  double w_x_z, w_x_eta, sol_num, sol_denom;
  double beta_new;

  for (int m = 1; m < num_lambda + 1; ++m) {
    if (m == num_lambda)
      lambda_seq[m] = 0;    // last solution should be unregularized
    else
      lambda_seq[m] = lambda_seq[m - 1] * eps_ratio;
    //std::cout << "\nm " << m << " lambda " << lambda_seq[m] << ", scaled: " << lambda_seq[m] * std_y << """ \n";

    /* Initialize the solution at this lambda using previous lambda solution */
    // We need to explicitly do this because we need to store all the solutions separately
    for (ull i = 0; i < n_params; ++i) {
      beta[i + n_params] = beta[i];
    }
    beta = beta + n_params;   // go to the next column
    for (ull i = 0; i < n_params; ++i) {
    }

    // log_scale += scale_update; scale = exp(log_scale);
    // std::cout  << "\n\n\n m: " << m << "\n";
    // std::cout << "-  log scale: " << log_scale << ", scale: " << scale << " update " << scale_update << "\n";

    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    n_iters[m] = 0;
    do {                                  // until Convergence of beta
      flag_beta_converged = 1;              // = 1 if beta converges
      //std::cout << "beta \n";

      /* iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_params; ++k) {

        sol_num = sol_denom = 0;          // TODO: You should optimize this so that we don't calculate the whole thing everytime
        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
          sol_num += (w[i] * X(i, k) * (z[i] - eta[i]));
          sol_denom += (w[i] * X(i, k) * X(i, k));
        }

        /* The intercept should not be regularized, and hence is calculated directly */
        if (intercept && k == 0) {
          beta_new = sol_num / sol_denom;   // TODO: Check

        } else {
          beta_new = soft_threshold(sol_num, lambda_seq[m] * alpha) /
                        (sol_denom + lambda_seq[m] * (1 - alpha));

        }
        // if any beta_k has not converged, we will come back for another cycle.
        if (fabs(beta_new - beta[k]) > threshold)
          flag_beta_converged = 0;

        beta[k] = beta_new;

        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
        }

      }   // end for: beta_k solution

      if (estimate_scale) {
        // calculate w and z again (beta & hence eta would have changed)
        compute_grad_response(NULL, NULL, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                              censoring_type, n_obs, transformed_dist, NULL);

        log_scale += scale_update; scale = exp(log_scale);
        // std::cout << "log scale: " << log_scale << ", scale: " << scale << " update " << scale_update << "\n";

        if (scale_update > threshold)
          flag_beta_converged = 0;
      }

      n_iters[m]++;
    } while ((n_iters[m] < max_iter) && (flag_beta_converged != 1));

    /* beta and eta will already contain their updated values since we calculate them in place */
    // calculate w and z again (beta & hence eta would have changed)
    compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                          censoring_type, n_obs, transformed_dist, NULL);
    // std::cout <<" w z - " << m << " \n";
    // for (ull i=0; i<n_obs; ++i) {
    //   std::cout << w[i] << " " << z[i] << "\n";
    // }
    // std::cout << "\n";
    // std::cout << "\n";

    out_scale[m] = scale;

  } // end for: lambda

  /* Compute the final log likelihood */
  loglik = compute_loglik(REAL(y), REAL(y) + n_obs, eta, scale,
                          censoring_type, n_obs, transformed_dist);

  /* Scale the coefs back to the original scale */
  if (standardize) {
    for (ull m = 0; m < num_lambda; ++m) {
      lambda_seq[m] = lambda_seq[m] * std_y;

      for (ull i = 0; i < n_vars; ++i) {
        out_beta(i, m) = out_beta(i, m) * std_y / std_x[i];
      }
    }
  }

  /* Provide values for the intercepts */
  if (standardize && intercept) {
    double *ptr = REAL(out_intercept), temp = 0;

    for (ull m = 0; m < num_lambda; ++m) {
      temp = 0.0;
      for (ull j = 0; j < n_vars; ++j) {
        temp += mean_x[j] * out_beta(j ,m);
      }

      ptr[m] = mean_y - temp;
    }
  }

  /* Free the temporary variables */
  delete [] eta;
  delete [] w;
  delete [] z;

  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta,
                            Rcpp::Named("lambda")       = out_lambda,
                            Rcpp::Named("n_iters")      = out_n_iters,
                            //Rcpp::Named("score")        = out_score,
                            Rcpp::Named("loglik")       = loglik,
                            Rcpp::Named("scale")        = out_scale,
                            Rcpp::Named("estimate_scale") = estimate_scale,
                            Rcpp::Named("error_status") = 0
                            );
}

  /* Convert the input to normalized form e_i = (log(y_i) - x_i' beta) / sigma */

static inline void get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *censoring_type)
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
IREG_DIST get_ireg_dist (Rcpp::String dist_str)
{
  if (strcmp("gaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_GAUSSIAN;
  if (strcmp("logistic", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOGISTIC;
  if (strcmp("extreme_value", dist_str.get_cstring()) == 0)
    return IREG_DIST_EXTREME_VALUE;
  if (strcmp("exponential", dist_str.get_cstring()) == 0)
    return IREG_DIST_EXPONENTIAL;
  if (strcmp("loggaussian", dist_str.get_cstring()) == 0)
    return IREG_DIST_LOG_GAUSSIAN;
  return IREG_DIST_UNKNOWN;
}

static inline double soft_threshold(double x, double lambda)
{
  double temp = fabs(x) - lambda;
  return (temp > 0)? ((x > 0)? 1: -1) * temp: 0;
}

static void standardize_x_y(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                            double *mean_x, double *std_x, double &mean_y,
                            double &std_y, bool intercept=false)
{
  double temp;
  ull count_y = 0, n_rows_y = y.size() / 2;

  /* Standardize y: mean and variance normalization */
  // Find mean if intercept needs to be fit, else set it to 0
  mean_y = std_y = 0;
  //for (ull i = 0; i < y.size(); ++i) {
  for (ull i = 0; i < n_rows_y; ++i) {
    if (y[i] == Rcpp::NA) continue;

    // only count if not NA
    count_y++;
    mean_y += y[i];
  }
  mean_y = mean_y / count_y;

  for (ull i = 0; i < n_rows_y; ++i) {
    if (y[i] == Rcpp::NA) continue;

    temp = y[i] - mean_y;
    std_y += temp * temp;
  }

  std_y = sqrt(std_y / count_y);
  // std_y = sqrt(std_y - mean_y * mean_y);

  // FIXME: why divide by sqrt(N)?
  for (ull i = 0; i < 2 * n_rows_y; ++i) {
    temp = (std_y * sqrt(count_y));      // FIXME: !!!
    if (intercept) {
      y[i] = (y[i] - mean_y) / temp;
    } else {
      y[i] = y[i] / temp;
    }
  }

  /* Mean and var normalize columns of X matrix */
  for (ull i = 0; i < X.ncol(); ++i) {
    mean_x[i] = std_x[i] = 0;
    for (ull j = 0; j < X.nrow(); ++j) {
      mean_x[i] += X(j, i);
    }
    mean_x[i] /= X.nrow();

    for (ull j = 0; j < X.nrow(); ++j) {
      temp = X(j, i) - mean_x[i];
      std_x[i] += temp * temp;
    }

    std_x[i] = sqrt(std_x[i] / X.nrow());
    // std_x[i] = sqrt(std_x[i] - mean_x[i] * mean_x[i]);

    for (ull j = 0; j < X.nrow(); ++j) {
      temp = (std_x[i] * sqrt(X.nrow()));
      if (intercept) {
        X(j, i) = (X(j, i) - mean_x[i]) / temp;         // FIXME: AS IN GLMNET!?
      } else {
        X(j, i) = X(j, i) / temp;         // FIXME: AS IN GLMNET!?
      }
      // For no intercept, return mean_x = 0;
      //mean_x[i] = 0;
    }
  }
}
