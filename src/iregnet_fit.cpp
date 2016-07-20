#include "iregnet.h"

#define BIG 1e15

void get_y_means (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *ym);
static inline double soft_threshold (double x, double lambda);
static void standardize_x (Rcpp::NumericMatrix X,
                           double *mean_x, double *std_x,
                           bool intercept);
static double get_init_var (double *ym, IREG_CENSORING *status, ull n, IREG_DIST dist);
static double get_init_intercept (double *mu, double *w, double *ym, ull n_obs);

double identity (double y)
{
  return y;
}

static double max(double a, double b)
{
  if(a > b)
    return a;
  return b;
}

static inline double
compute_lambda_max(Rcpp::NumericMatrix X, double *w, double *z, bool intercept,
									 double &alpha, ull n_vars, ull n_obs);

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
 * This follows the math as outlined
 */
// [[Rcpp::export]]
Rcpp::List fit_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y,
                   Rcpp::String family,   Rcpp::NumericVector lambda_path,
                   Rcpp::IntegerVector out_status,
                   bool intercept,        double alpha,
                   double scale_init,     bool estimate_scale,
                   bool unreg_sol,        bool flag_standardize_x,
                   double max_iter,       double threshold,
                   int num_lambda,        double eps_lambda)
{
  /* Initialise some helper variables */
  ull n_obs, n_vars;
  IREG_DIST orig_dist, transformed_dist;  // Orig dist is the one that is initially given
                                          // transformed_dist will be the dist of transformed output variables.
																					// Eg- orig_dist = "loglogistic", transformed_dist="logistic", with transform_y=log
  double (*transform_y) (double y);
	double scale;

  n_obs  = X.nrow();
  n_vars = X.ncol();  // n_vars is the number of variables corresponding to the coeffs of X
  orig_dist = get_ireg_dist(family);

  /* TODO: Validate all the arguments again */
  // if ((alpha > 1 || alpha < 0) || (IREG_DIST_UNKNOWN == orig_dist) ||
  //     (y.nrow() != n_obs)) {
  //   return Rcpp::List::create(Rcpp::Named("error_status") = -1);
  // }


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
      estimate_scale = false;
      break;
  }

  for (ull i = 0; i < 2 * n_obs; ++i) {     // Rcpp::NumericMatrix is unrolled col major
    y[i] = transform_y(y[i]);
  }

  /*
   * Standardize functions:
   */
  double *mean_x = new double [n_vars];
  double *std_x = new double [n_vars];

  /* NOTE:
   * we have scaled y and x, so you need to scale the obtained lambda values,
   * and coef (beta) values back to the original scale before returning them.
   */
  if (flag_standardize_x) {
    standardize_x(X, mean_x, std_x, intercept);      // FIXME: so that values are always estimated as for intercepts
  }

  /* Create output variables */
	if (lambda_path.size() > 0) {
		num_lambda = lambda_path.size();
	}
  Rcpp::NumericMatrix out_beta(n_vars, num_lambda);       // will contain the entire series of solutions
  Rcpp::IntegerVector out_n_iters(num_lambda);
  Rcpp::NumericVector out_lambda(num_lambda);
  Rcpp::NumericVector out_scale(num_lambda);
  Rcpp::NumericVector out_loglik(num_lambda);

	/* use given values for the lambda path */
	if (lambda_path.size() > 0) {
		for(ull i = 0; i < num_lambda; ++i)
			out_lambda[i] = lambda_path[i];
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
	double *ym = new double [n_obs];		// the mean values for y's
  IREG_CENSORING *status;


  /* get censoring types of the observations */
  if (out_status.size() == 0) {
    status = new IREG_CENSORING [n_obs];
    get_censoring_types(y, status); // NANs denote censored observations in y, NOTE: y will be modified!
  }
  else {
    status = (IREG_CENSORING *) &out_status[0];
  }
  get_y_means(y, status, ym);

	/* SCALE RULES:
	 * Whether or not it is estimated depends on estimate_scale. For exponential, this is forced to False and scale fixed to 1. // TODO
	 *
	 * If you provide no scale, a starting value will be calculated
	 * If you provide a scale, it will be used as the initial value
	 */
  if (scale_init == Rcpp::NA) {
    scale = get_init_var(ym, status, n_obs, transformed_dist);
		// scale = sqrt(scale);		// use 2 * sqrt(var) as in survival
		scale = 2 * sqrt(scale);		// use 2 * sqrt(var) as in survival
		// TODO: find corner cases where you need to take 2 * sqrt(var) as in survival
	} else {
		scale = scale_init;
	}
	log_scale = log(scale);

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

  for (int m = 0; m < num_lambda; ++m) {
		/* Compute the lambda path */
		if (lambda_path.size() == 0) {

			/* Do an initial fit with lambda set to BIG, will fit scale and intercept if applicable */
			if (m == 0) {
				lambda_max_unscaled = lambda_seq[0] = 1e35;
			}

			/* Calculate lambda_max using intial scale fit */
			if (m == 1) {
				lambda_seq[m] = compute_lambda_max(X, w, z, intercept, alpha, n_vars, n_obs);
				lambda_max_unscaled = lambda_seq[m] * scale * scale;
				// lambda_max_unscaled = lambda_seq[m] = lambda_seq[m] * scale * scale;

			/* Last solution should be unregularized if the flag is set */
			} else if (m == num_lambda - 1 && unreg_sol == true)
    	  lambda_seq[m] = 0;

			/* All other lambda calculated */
      else if (m > 1) {
        // lambda_seq[m] = lambda_seq[m - 1] * eps_ratio;
        lambda_seq[m] = lambda_max_unscaled * pow(eps_ratio, m-1) / scale / scale;
      }
		}

    /* Initialize the solution at this lambda using previous lambda solution */
    // We need to explicitly do this because we need to store all the solutions separately
		if (m != 0) {
			for (ull i = 0; i < n_vars; ++i) {
				beta[i + n_vars] = beta[i];
			}
			beta = beta + n_vars;   // go to the next column
		}

    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    n_iters[m] = 0;
    do {                                  // until Convergence of beta

      flag_beta_converged = 1;            // = 1 if beta converges
			old_scale = scale;

      // IRLS: Reweighting step: calculate w and z again (beta & hence eta would have changed)  TODO: make dg ddg, local so that we can save computations?
      compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                            status, n_obs, transformed_dist, NULL);

      /* iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_vars; ++k) {
        sol_num = sol_denom = 0;
        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
          sol_num += (w[i] * X(i, k) * (z[i] - eta[i])) / n_obs;
          sol_denom += (w[i] * X(i, k) * X(i, k)) / n_obs;
        }

        // Note: The signs given in the coxnet paper are incorrect, since the subdifferential should have a negative sign.
        sol_num *= -1; sol_denom *= -1;

        /* The intercept should not be regularized, and hence is calculated directly */
        if (intercept && k == 0) {
          beta_new = sol_num / sol_denom;

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
        log_scale += scale_update; scale = exp(log_scale);
				// scale the lambda value according to current scale unless you are at unregularized sol
				// done to match values with Glment for Gaussian, no censoring
        if ((m != num_lambda - 1 || unreg_sol == false) && lambda_path.size() == 0 && m > 1)
          lambda_seq[m] = lambda_max_unscaled * pow(eps_ratio, m-1) / scale / scale;    // FIXME: Scale! :O

        // if (fabs(scale - old_scale) > threshold) {		// TODO: Maybe should be different for sigma?
        if (fabs(scale - old_scale) > 1e-4) {		// TODO: Maybe should be different for sigma?
          flag_beta_converged = 0;
        }
      }

      n_iters[m]++;
    } while ((flag_beta_converged != 1) && (n_iters[m] < max_iter));

    out_scale[m] = scale;
		if (m == 0)
			scale_init = scale;
  } // end for: lambda

  /* Scale the coefs back to the original scale */
  if (flag_standardize_x) {
    for (ull m = 0; m < num_lambda; ++m) {

      for (ull i = int(intercept); i < n_vars; ++i) {	// no standardization for intercept
        out_beta(i, m) = out_beta(i, m) / std_x[i];
      }
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

  return Rcpp::List::create(Rcpp::Named("beta")         = out_beta,
                            Rcpp::Named("lambda")       = out_lambda,
                            Rcpp::Named("num_lambda")   = num_lambda,
                            Rcpp::Named("n_iters")      = out_n_iters,
                            Rcpp::Named("loglik")       = out_loglik,
                            Rcpp::Named("scale")        = out_scale,
                            Rcpp::Named("estimate_scale") = estimate_scale,
                            Rcpp::Named("scale_init")   = scale_init,
                            Rcpp::Named("error_status") = 0
                            );
}

void get_y_means (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *ym)
{
  if (!ym || !status)
    return;

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
      // std::cout<<"status, ym: " << status[i] << " " << ym[i] << "\n";
  } // end for
}

void get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *status)
{
  for (ull i = 0; i < y.nrow(); ++i) {
    if (y(i, 0) == Rcpp::NA) {

      if (y(i, 1) == Rcpp::NA)
        status[i] = IREG_CENSOR_INVALID;    // invalid data
      else {
        status[i] = IREG_CENSOR_LEFT;       // left censoring
        y(i, 0) = y(i, 1); // NOTE: We are putting the value in the left col, survival style!
      }
      continue;
    }

    if (y(i, 1) == Rcpp::NA)                // right censoring
      status[i] = IREG_CENSOR_RIGHT;
    else if (y(i, 0) == y(i, 1))
      status[i] = IREG_CENSOR_NONE;         // no censoring
		else
      status[i] = IREG_CENSOR_INTERVAL;     // interval censoring

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

static inline double soft_threshold (double x, double lambda)
{
  double temp = fabs(x) - lambda;
  return (temp > 0)? ((x > 0)? 1: -1) * temp: 0;
}

static void standardize_x (Rcpp::NumericMatrix X,
                           double *mean_x, double *std_x,
                           bool intercept)
{
  double temp;

  for (ull i = int(intercept); i < X.ncol(); ++i) {	// don't standardize intercept col.
    mean_x[i] = std_x[i] = 0;
    for (ull j = 0; j < X.nrow(); ++j) {
      mean_x[i] += X(j, i);
    }
    mean_x[i] /= X.nrow();

    for (ull j = 0; j < X.nrow(); ++j) {
      temp = X(j, i) - mean_x[i];
      std_x[i] += temp * temp;
    }

    // not using (N-1) in denominator to agree with glmnet
    std_x[i] = sqrt(std_x[i] / X.nrow());

    for (ull j = 0; j < X.nrow(); ++j) {
      X(j, i) = X(j, i) / std_x[i];
    }
  }
}

static double get_init_var (double *ym, IREG_CENSORING *status, ull n, IREG_DIST dist)
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

		default:
		break;
	}

	return var;
}

static double get_init_intercept (double *mu, double *w, double *ym, ull n_obs)
{
	double intercept = 0;

	for (ull i = 0; i < n_obs; ++i) {
		intercept += (mu[i] - w[i] * ym[i]);
	}
	intercept /= n_obs;

	return intercept;
}

static inline double
compute_lambda_max(Rcpp::NumericMatrix X, double *w, double *z, bool intercept,
									 double &alpha, ull n_vars, ull n_obs)
{
	double lambda_max = -1; // calculated values are always non-negative
	for (ull j = int(intercept); j < n_vars; ++j) {   // dont include intercept col in lambda calc.
	  double temp = 0;

	  for (ull i = 0; i < n_obs; ++i) {
	    temp += (w[i] * X(i, j) * z[i]);
	  }
		temp = fabs(temp);
	  lambda_max = (lambda_max > temp)? lambda_max: temp;
	}
	lambda_max /= (n_obs * max(alpha, 1e-3));  // prevent divide by zero

	return lambda_max;
}
