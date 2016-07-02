#include "iregnet.h"

#define BIG 1e35

static inline void get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *yy);
static inline double soft_threshold (double x, double lambda);
static void standardize_x (Rcpp::NumericMatrix X,
                           double *mean_x, double *std_x,
                           bool intercept);
static double get_init_var (double *yy, IREG_CENSORING *status, ull n, IREG_DIST dist);
static double get_init_intercept (double *mu, double *w, double *yy, ull n_obs);
static inline void fit_scale_intercept(double *w, double *z, double &scale_update, Rcpp::NumericMatrix X, Rcpp::NumericMatrix y, ull &n_obs, double *eta,
																		 	 double &scale, IREG_CENSORING *status, IREG_DIST &transformed_dist, double *mu, double *beta, double &threshold,
																			 double &max_iter, double &log_scale, bool &intercept);
double identity (double y)
{
  return y;
}

/*
 * TODO: Better initial value for mean using glim trick.. Look at survival!
 * TODO: As a result (?) you should have a good value for the lambda path
 */

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
                   Rcpp::NumericVector lambda_path, bool intercept,
                   double scale,      bool estimate_scale = true,
                   bool flag_standardize_x = false,
                   double max_iter = 1000,  double threshold = 1e-4,
                   int num_lambda = 100,  double eps_lambda = 0.0001,   // TODO: depends on nvars vs nobs
                   int flag_debug = 0)
{
  /* Initialise some helper variables */
  ull n_obs, n_vars, n_cols_x, n_cols_y, n_params;
  IREG_DIST orig_dist, transformed_dist;  // Orig dist is the one that is initially given
                                          // transformed_dist will be the dist of transformed output variables.
																					// Eg- orig_dist = "loglogistic", transformed_dist="logistic", with transform_y=log
  double (*transform_y) (double y);

  n_obs  = X.nrow();
  n_cols_x = X.ncol();
  n_cols_y = y.ncol();
  orig_dist = get_ireg_dist(family);

  n_vars = X.ncol();  // n_vars is the number of variables corresponding to the coeffs of X
  n_params = n_vars; // n_params is the number of parameters

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
      estimate_scale = false;
      break;
  }

  // TODO: unnecessary computations
  for (ull i = 0; i < 2 * n_obs; ++i) {     // Rcpp::NumericMatrix is unrolled col major
    y[i] = transform_y(y[i]);
  }

  /*
   * Standardize functions:
   */
  double *mean_x = new double [n_vars], mean_y;
  double *std_x = new double [n_vars], std_y;

  /* NOTE:
   * we have scaled y and x, so you need to scale the obtained lambda values,
   * and coef (beta) values back to the original scale before returning them.
   */
  if (flag_standardize_x) {
    standardize_x(X, mean_x, std_x, intercept);      // FIXME: so that values are always estimated as for intercepts
    // standardize_x_y(X, y, mean_x, std_x, mean_y, std_y, true);
    // std::cout << "y\n" << y << "\nx\n" << X;
    // std::cout << "mean_y " << mean_y << "std_y " << std_y << "mean_x:\n";
    // for (ull i=0; i<n_vars; ++i) {
    //   std::cout << i << " " << mean_x[i] << " " << std_x[i] << "\n";
    // }
  }

  /* Create output variables */
	if (lambda_path.size() > 0) {
		num_lambda = lambda_path.size();
		// std::cout<<"Using given values of lambda\n";
	}
  Rcpp::NumericMatrix out_beta(n_params, num_lambda + 1);       // will contain the entire series of solutions
  Rcpp::IntegerVector out_n_iters(num_lambda + 1);
  Rcpp::NumericVector out_lambda(num_lambda + 1);
  Rcpp::NumericVector out_scale(num_lambda + 1);
  Rcpp::NumericVector out_loglik(num_lambda + 1);

	/* use given values for the lambda path */
	if (lambda_path.size() > 0) {
		for(ull i = 0; i < num_lambda; ++i)
			out_lambda[i] = lambda_path[i];
		out_lambda[num_lambda] = 0;
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
	double *mu = new double [n_obs];
	double *yy = new double [n_obs];		// the values for y's

	double scale_init = 0;

  /* get censoring types of the observations */ /* TODO: Incorporate survival style censoring */
  IREG_CENSORING status[n_obs];
  get_censoring_types(y, status, yy); // NANs denote censored observations in y

	/* SCALE RULES:
	 * Whether or not it is estimated depends on estimate_scale. For exponential, this is forced to False and scale fixed to 1. // TODO
	 *
	 * If you provide no scale, a starting value will be calculated
	 * If you provide a scale, it will be used as the initial value
	 */
  if (scale == Rcpp::NA) {
    scale = get_init_var(yy, status, n_obs, transformed_dist);
		// scale = sqrt(scale);		// use 2 * sqrt(var) as in survival
		scale = 2 * sqrt(scale);		// use 2 * sqrt(var) as in survival
		// TODO: find corner cases where you need to take 2 * sqrt(var) as in survival
	}
	log_scale = log(scale);
	// std::cout<< "Init scale value is : " << scale << " logscale: " << log_scale << "\n";

  /* Initalize the pathwise solution
   * We will always start with lambda_max, at which beta = 0, eta = 0.
   * If some value of lambda is supplied, we will stop iterations at that value.
   */

  // set beta = 0 and eta = 0
	beta = REAL(out_beta);
  for (ull i = 0; i < n_params; ++i) {
    // sets only the first solution (first col of out_beta) to 0
    beta[i] = 0;
  }
  n_iters[0] = 0; out_scale[0] = scale;

  for (ull i = 0; i < n_obs; ++i) {
    eta[i] = w[i] = z[i] = 0;
  }


  /* ************************************************************************************ */
	/***** fit0: Fit only intercept and scale **************************************************************************/
	/* 	 * except when scale doesn't need to be estimated, or initial value already given
	 * TODO: Mean-only should only fit scale? TODO: What if we have no intercept -> fit only scale.?
	 */
	if (estimate_scale && intercept) {
		fit_scale_intercept(w, z, scale_update, X, y, n_obs, eta, scale, status, transformed_dist, mu, beta, threshold, max_iter, log_scale, intercept);
	}
	scale_init = scale;

	/***** END fit0 ************************************************************************************************/

  // Calculate w and z right here!
  compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,
                        status, n_obs, transformed_dist, mu);

  // Calculate lambda_max
  // First, start with lambda_max = BIG (really really big), so that eta and beta are surely 0
	if (lambda_path.size() == 0) {
  	lambda_seq[0] = -BIG;
  	for (ull j = int(intercept); j < n_vars; ++j) {   // dont include intercept col in lambda calc.
  	  double temp = 0;

  	  for (ull i = 0; i < n_obs; ++i) {
  	    temp += (w[i] * X(i, j) * z[i]);
  	  }
  	  //temp = fabs(temp) / (n_obs * alpha);
			temp = fabs(temp);
  	  lambda_seq[0] = (lambda_seq[0] > temp)? lambda_seq[0]: temp;
  	}
		lambda_seq[0] /= (n_obs * alpha);
	}

	// beta[0] = get_init_intercept(mu, w, yy, n_obs); // set initial value of intercept using glim trick (survival)
	// std::cout << beta[0] << " <<--------- beta[0]\n";

  // TODO: you should only set lambda_max = Inf, and let it calc beta = eta = 0 itself.

	/******************************************************************************************/
  /* Iterate over grid of lambda values */
  double eps_ratio = std::pow(eps_lambda, 1.0 / (num_lambda-1));
  bool flag_beta_converged = 0;
  double sol_num, sol_denom;
  double beta_new;
  double old_scale;
	double INDEX = 2000;

  // for (int m = 0; m < 2; ++m) {
  for (int m = 0; m < num_lambda + 1; ++m) {
		// TODO: Fix lambda_max and verify that rest of values are correct
		// TODO: The eps_ratio depends on some conditions, code them
		if (lambda_path.size() == 0) {
    	if (m == num_lambda)
    	  lambda_seq[m] = 0;    // last solution should be unregularized
    	else if (m != 0)
    	  lambda_seq[m] = lambda_seq[m - 1] * eps_ratio;
		}

    /* Initialize the solution at this lambda using previous lambda solution */
    // We need to explicitly do this because we need to store all the solutions separately
		if (m != 0) {
			for (ull i = 0; i < n_params; ++i) {
				beta[i + n_params] = beta[i];
			}
			beta = beta + n_params;   // go to the next column
		}

		old_scale = scale;

    /* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    n_iters[m] = 0;
    do {                                  // until Convergence of beta
			if (m == INDEX) {
				 for (int i=0; i<n_obs; ++i) {
				 	std::cout << "w: " << w[i] << " , z: " << z[i] << "\n";
				 }
				 Rcpp::NumericVector aa = out_beta(Rcpp::_, m);
				 std::cout << "lambda: " << lambda_seq[m] << "\n";
				 std::cout << "out_beta starting: " << aa << "\n";
			}

      flag_beta_converged = 1;              // = 1 if beta converges

      /* iterate over beta elementwise and update using soft thresholding solution */
      for (ull k = 0; k < n_params; ++k) {
        sol_num = sol_denom = 0;
        for (ull i = 0; i < n_obs; ++i) {
				  // if (k != 0)
					// 	eta[i] = eta[i] - X(i, k) * beta[0];
          eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
          sol_num += (w[i] * X(i, k) * (z[i] - eta[i])) / n_obs;
          sol_denom += (w[i] * X(i, k) * X(i, k)) / n_obs;
					if (m == INDEX)
						std::cout << eta[i] << ", ";
        }
				if (m == INDEX)
					std::cout << "\n";

        /* The intercept should not be regularized, and hence is calculated directly */
        if (intercept && k == 0) {
          beta_new = sol_num / sol_denom;

        } else {
          beta_new = soft_threshold(sol_num, lambda_seq[m] * alpha) /
                        (sol_denom + lambda_seq[m] * (1 - alpha));
					if (m == INDEX)
						std::cout << "sols: " << sol_num << " " << sol_denom << " " << beta_new << "\n \n";

        }
        // if any beta_k has not converged, we will come back for another cycle.
        if (fabs(beta_new - beta[k]) > threshold)
          flag_beta_converged = 0;

        beta[k] = beta_new;

        for (ull i = 0; i < n_obs; ++i) {
          eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
				  // if (k != 0)
					// 	eta[i] = eta[i] + X(i, k) * beta[0];
        }

      }   // end for: beta_k solution

      if (estimate_scale) {
        // std::cout << "log scale: " << log_scale << ", scale: " << scale << " update " << scale_update << "\n";
        log_scale += scale_update; scale = exp(log_scale);

        if (fabs(scale - old_scale) > threshold) {		// TODO: Maybe should be different for sigma?
          flag_beta_converged = 0;
          // calculate w and z again (beta & hence eta would have changed)  TODO: make dg ddg, local so that we can save computations?
        }
				old_scale = scale;
      }
      compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                            status, n_obs, transformed_dist, NULL);

      n_iters[m]++;
			if (m == INDEX) {
				 Rcpp::NumericVector aa = out_beta(Rcpp::_, m);
				 std::cout << "out_beta ending: " << aa << "\n";
			}
    } while ((n_iters[m] < max_iter) && (flag_beta_converged != 1));

    /* beta and eta will already contain their updated values since we calculate them in place */
    // calculate w and z again (beta & hence eta would have changed)
    out_loglik[m] = compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
                          status, n_obs, transformed_dist, NULL);

    out_scale[m] = scale;

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
  delete [] w;
  delete [] z;

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

static inline void get_censoring_types (Rcpp::NumericMatrix &y, IREG_CENSORING *status, double *yy)
{
  for (ull i = 0; i < y.nrow(); ++i) {
    if (y(i, 0) == Rcpp::NA) {

      if (y(i, 1) == Rcpp::NA)
        status[i] = IREG_CENSOR_INVALID;    // invalid data
      else {
        status[i] = IREG_CENSOR_LEFT;       // left censoring
				yy[i] = y(i, 1);
			}

    } else {
      if (y(i, 1) == Rcpp::NA) {                      // right censoring
        status[i] = IREG_CENSOR_RIGHT;
				yy[i] = y(i, 0);
			}

      else if (y(i, 0) == y(i, 1)) {
        status[i] = IREG_CENSOR_NONE;     // no censoring
				yy[i] = y(i, 1);
			} else {
        status[i] = IREG_CENSOR_INTERVAL; // interval censoring
				yy[i] = (y(i, 1) + y(i, 0)) / 2;
			}
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

    std_x[i] = sqrt(std_x[i] / X.nrow());
    // std_x[i] = sqrt(std_x[i] - mean_x[i] * mean_x[i]);

    for (ull j = 0; j < X.nrow(); ++j) {
      temp = std_x[i];
      //temp = (std_x[i] * sqrt(X.nrow()));
      // X(j, i) = X(j, i) / temp;         // FIXME: AS IN GLMNET!?
    }
  }
}

static double get_init_var (double *yy, IREG_CENSORING *status, ull n, IREG_DIST dist)
{
	double mean = 0, var = 0;

	for (int i = 0; i < n; ++i) {
		mean += yy[i];
	}
	mean = mean / n;

	for (int i = 0; i < n; ++i) {
		var += pow((yy[i] - mean), 2);
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

static double get_init_intercept (double *mu, double *w, double *yy, ull n_obs)
{
	double intercept = 0;

	for (ull i = 0; i < n_obs; ++i) {
		intercept += (mu[i] - w[i] * yy[i]);
	}
	intercept /= n_obs;

	return intercept;
}

static inline void fit_scale_intercept(double *w, double *z, double &scale_update, Rcpp::NumericMatrix X, Rcpp::NumericMatrix y, ull &n_obs, double *eta,
																		 	 double &scale, IREG_CENSORING *status, IREG_DIST &transformed_dist, double *mu, double *beta, double &threshold,
																			 double &max_iter, double &log_scale, bool &intercept)
{
			double n_iters_fit0 = 0;
  		bool flag_beta_converged = 0;
  		double sol_num, sol_denom;
  		double beta_new;
  		double old_scale;

  		compute_grad_response(w, z, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,
  		                      status, n_obs, transformed_dist, mu);

			old_scale = scale;

    	/* CYCLIC COORDINATE DESCENT: Repeat until convergence of beta */
    	n_iters_fit0 = 0;
    	do {                                  // until Convergence of beta
    	  flag_beta_converged = 1;              // = 1 if beta converges

    	  /* iterate over beta elementwise and update using soft thresholding solution */
    	  //for (ull k = 0; k < n_params; ++k) {
				ull k = 0;		// only fit the intercept

    	    sol_num = sol_denom = 0;          // TODO: You should optimize this so that we don't calculate the whole thing everytime
    	    for (ull i = 0; i < n_obs; ++i) {
    	      eta[i] = eta[i] - X(i, k) * beta[k];  // calculate eta_i without the beta_k contribution
    	      sol_num += (w[i] * X(i, k) * (z[i] - eta[i])) / n_obs;
    	      sol_denom += (w[i] * X(i, k) * X(i, k)) / n_obs;
    	    }

    	    /* The intercept should not be regularized, and hence is calculated directly */
    	    if (intercept && k == 0)
    	      beta_new = sol_num / sol_denom;

    	    // if any beta_k has not converged, we will come back for another cycle.
    	    if (fabs(beta_new - beta[k]) > threshold)
    	      flag_beta_converged = 0;

    	    beta[k] = beta_new;

    	    for (ull i = 0; i < n_obs; ++i) {
    	      eta[i] = eta[i] + X(i, k) * beta[k];  // this will contain the new beta_k
    	    }

    	  // }   // end for: beta_k solution

    	  //if (estimate_scale) {
    	    log_scale += scale_update; scale = exp(log_scale);

    	      compute_grad_response(NULL, NULL, &scale_update, REAL(y), REAL(y) + n_obs, eta, scale,     // TODO:store a ptr to y?
    	                            status, n_obs, transformed_dist, NULL);
    	    if (fabs(scale - old_scale) > threshold) {		// TODO: Maybe should be different for sigma?
    	      flag_beta_converged = 0;
    	      // calculate w and z again (beta & hence eta would have changed)  TODO: make dg ddg, local so that we can save computations?
    	    }
					old_scale = scale;
    	  //}

    	  n_iters_fit0++;
    	} while ((n_iters_fit0 < max_iter) && (flag_beta_converged != 1));
    // std::cout << "fit0: log scale: " << log_scale << ", scale: " << scale << " beta[0] " << beta[0]  << " n_iters " << n_iters_fit0 << "\n";
}
