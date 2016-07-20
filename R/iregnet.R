# iregnet funciton
# USAGE
# Authors
# Borrowed parts from survival survreg and glmnet

# TODO: optimization params need to be implemented
# TODO: checks like - can't provide scale with exp
iregnet <- function(x, y,
                    family=c("gaussian", "logistic", "loggaussian", "extreme_value", "exponential"), alpha=1,
                    lambda=double(0), num_lambda=100, flag_debug=0, intercept=T, standardize=F, scale_init=NA, estimate_scale=T,
										maxiter=1000, threshold=1e-4, unreg_sol=T, eps_lambda=NA) {

  # Parameter validation ===============================================
  # alpha should be between 0 and 1
  if (0 > alpha || alpha > 1) {
    stop("alpha should be between 0 and 1")
  }

	if (estimate_scale == F && is.na(scale_init))
		stop("Value of scale required if scale is not estimated")

  # family should be one of those listed
  family <- match.arg(family)

  # x should be a matrix with atleast two columns (TODO: why glmnet requires 2 or more cols?)
    # TODO: NAs in x? should be numeric
  if (!is.matrix(x)) stop("x should be a matrix with 2 or more columns")
  if (ncol(x) <= 1) stop("x should be a matrix with 2 or more columns")

  n_obs  <- nrow(x)   # LIST
  n_vars <- ncol(x)   # LIST

  # y should be a matrix with 2 columns correspoding to the left and right times
  # NA or Inf/-Inf can be used for censored entries
  if (is.matrix(y)){

  # OR y can be a Surv object from the survival library
  } else if (survival::is.Surv(y)) {
  }
  dim_y <- dim(y);      # LIST
  if (is.null(dim_y) | dim_y[2] != 2)
   stop("y should be a matrix with 2 columns");

  # length of x and y should be the same
  if (dim_y[1] != n_obs)
    stop("number of rows in x and y don't match");

	# Append col of 1's for the intercept
	if (intercept)
		x <- cbind(rep(1, n_obs), x)

  if (is.na(eps_lambda))
    eps_lambda <- ifelse(n_obs < n_vars, 0.01, 0.0001)
  # Call the actual fit method
  fit_cpp(x, y, family, alpha, lambda_path=lambda, num_lambda=num_lambda, intercept=intercept, scale_init=scale_init, max_iter=maxiter,
            flag_standardize_x = standardize, threshold=threshold, estimate_scale=estimate_scale, unreg_sol=unreg_sol, eps_lambda=eps_lambda);
}
