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
  stopifnot_error("alpha should be between 0 and 1", 0 <= alpha, alpha <= 1)

	if (estimate_scale == F && is.na(scale_init))
		stop("Value of scale required if scale is not estimated")

  # family should be one of those listed
  family <- match.arg(family)

  # TODO: NAs in x? should be numeric
  stopifnot_error("x should be a matrix with 2 or more columns", is.matrix(x), ncol(x) > 1)

  n_obs  <- nrow(x)
  n_vars <- ncol(x)

  # y should be a matrix with 2 columns correspoding to the left and right times
  # NA or Inf/-Inf can be used for censored entries
  stopifnot_error("y should be a 2 column matrix, or a Surv object", is.matrix(y) || survival::is.Surv(y))

  if (is.matrix(y)){
    stopifnot_error("y should be a 2 column matrix with nrow(y) = nrow(x)", ncol(y) == 2, nrow(y) == n_obs)

  } else if (survival::is.Surv(y)) {
  }

	# Append col of 1's for the intercept
	if (intercept)
		x <- cbind(rep(1, n_obs), x)

  if (is.na(eps_lambda))
    eps_lambda <- ifelse(n_obs < n_vars, 0.01, 0.0001)
  stopifnot_error("eps_lambda should be in [0, 1)", 0 <= eps_lambda && eps_lambda < 1)

  # Call the actual fit method
  fit <- fit_cpp(x, y, family, alpha, lambda_path=lambda, num_lambda=num_lambda, intercept=intercept, scale_init=scale_init, max_iter=maxiter,
                 flag_standardize_x = standardize, threshold=threshold, estimate_scale=estimate_scale, unreg_sol=unreg_sol, eps_lambda=eps_lambda);

  fit$call <- match.call()
  fit
}

# like stopifnot, but with a custom error message
stopifnot_error <- function(err_message, ...)
{
  n <- length(ll <- list(...))
  for (i in 1:n)
    if (!(is.logical(r <- ll[[i]]) && !anyNA(r) && all(r))) {
      stop(err_message)
    }
}
