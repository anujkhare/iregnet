# iregnet funciton
# USAGE
# Authors
# Borrowed parts from survival survreg and glmnet

iregnet <- function(x, y,
                    family=c("gaussian", "logistic", "loggaussian", "loglogistic", "extreme_value", "exponential"), alpha=1,
                    lambda=double(0), num_lambda=100, intercept=T, standardize=F, scale_init=NA, estimate_scale=T,
                    maxiter=1000, threshold=1e-4, unreg_sol=T, eps_lambda=NA, debug=F) {

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

  status <- integer(0) # used for censoring, if y is matrix, calculated in C++ code
  if (survival::is.Surv(y)) {
    status <- get_status_from_surv(y)
    y <- as.matrix(y[, 1:(ncol(y)-1)])
  } else {
    stopifnot_error("y should be a 2 column matrix with nrow(y) = nrow(x)", ncol(y) == 2, nrow(y) == n_obs)
  }
  stopifnot_error("y should be positive for the given family",
                  !(family %in% c('loglogistic', 'loggaussian', 'weibull') && any(y[!is.na(y)]<0)))

  # Append col of 1's for the intercept
  if (intercept)
    x <- cbind(rep(1, n_obs), x)

  if (is.na(eps_lambda))
    eps_lambda <- ifelse(n_obs < n_vars, 0.01, 0.0001)
  stopifnot_error("eps_lambda should be in [0, 1)", 0 <= eps_lambda && eps_lambda < 1)

  # Call the actual fit method
  fit <- fit_cpp(x, y, family, alpha, lambda_path=lambda, num_lambda=num_lambda, intercept=intercept,
                 out_status=status, scale_init=scale_init, max_iter=maxiter, threshold=threshold,
                 flag_standardize_x=standardize, estimate_scale=estimate_scale, unreg_sol=unreg_sol,
                 eps_lambda=eps_lambda, debug=debug);

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

get_status_from_surv <- function(s)
{
  type <- attr(s, 'type')

  stopifnot_error("Unsupported censoring type from Surv", type == 'left' || type == 'right' ||
                                                          type == 'interval' || type == 'interval2')
  # right censored: 0, none: 1, left: 2, interval: 3
  status <- s[, ncol(s)]
  if (type == 'left')
      status[status == 0] <- 2

  return(status)
}
