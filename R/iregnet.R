#' @title Fit interval censored AFT models with elastic net regularization
#'
#' @description
#' Fit accelerated failure time models using interval censored data via
#' elastic net penalized maximum likeihood. Solutions are computed using
#' coordinate descent for a path of values of the regularization parameter
#' lambda. Supports gaussian, logistic and extreme value distributions.
#'
#' @param x Input matrix of covariates with dimension n_obs * n_vars, with
#' \eqn{nvars \ge 2}. Sparse matrices are not supported.
#'
#' @param y Response variable. It can take two forms: \itemize{
#' \item 2 column real matrix with NAs denoting a censored observation
#' \item \code{\link{Surv}} object. Supported \code{types} of \code{Surv}
#' are 'left', 'right', 'interval' and 'interval2'.
#' }
#'
#' @param family The distribution to fit. It can be one of "gaussian", "logistic",
#' "loggaussian", "loglogistic", "extreme_value", "exponential". Partial matching
#' is allowed.
#' \cr \emph{Default: "gaussian"}
#'
#' @param alpha Elastic net mixing parameter, with \eqn{0 \le \alpha \le 1}. The
#' elastic net penalty is defined as in \code{glmnet}:
#' \deqn{0.5 * (1-\alpha) \|\beta\|_2^2 | + \alpha \|\beta\|_1}
#' alpha=1 is the lasso penalty, and alpha=0 is ridge penalty.
#' \cr \emph{Default: 1}
#'
#' @param lambda Vector containing the path of \strong{decreasing} regularization
#' parameter lambda values. \cr If not supplied, the function will calculate a
#' lambda path of length \code{num_lambda} itself. \strong{NOTE:} The lambda
#' values are scaled because of the nuisance parameter, and hence not directly
#' comparable to those of other packages like \code{glmnet}.
#' \cr \emph{Default: \code{NA}}
#'
#' @param num_lambda The number of lambda values calculated by the function.
#' Ignored if \code{lambda} is supplied by the user.
#' \cr \emph{Default: 100}
#'
#' @param intercept \code{TRUE} if an intercept is to be fit, otherwise
#' \code{FALSE}. Intercept is calculated by appending a column of 1s to \code{x}.
#' \cr \emph{Default: \code{TRUE}}
#'
#' @param standardize \code{TRUE} if \code{x} must be standardized before fit,
#' otherwise \code{FALSE}. Calculated \code{beta} are scaled to original scale
#' before returning.
#' \cr \emph{Default: \code{TRUE}}
#'
#' @param scale_init Initial value of the scale parameter to use. If not supplied,
#' a suitable value is calculated depending on the distribution.
#' \cr \emph{Default: NA}
#'
#' @param estimate_scale \code{TRUE} if \code{scale} is to be estimated. To
#' use a fixed value \code{scale0} for \code{scale}, set \code{scale_init=scale0
#' , estimate_scale=FALSE}. \emph{See examples.}
#' \cr \emph{Default: \code{TRUE}}
#'
#' @param maxiter Maximum number of iterations over data per lambda value.
#' \cr \emph{Default: 1e3}
#'
#' @param threshold The convergence threshold for coordinate descent. The inner
#' loop continues until the absolute update in each parameter is greater than
#' \code{threshold}.
#' \cr \emph{Default: 1e-4}
#'
#' @param eps_lambda The ratio of the minimum value of \code{lambda} to the
#' (calculated) maximum value, in case no lambda is supplied. \code{num_lambda}
#' \code{lambda} values are calculated between \code{lambda_max} and
#' \code{lambda_min} on the log scale.
#' \cr \emph{Default: 0.0001 if n_vars < n_obs, 0.1 otherwise.}
#'
#' @param unreg_sol \code{TRUE} if the final solution computed must be
#' unregularized. Overwritten to \code{FALSE} if n_vars > n_obs.
#' \cr \emph{Default: \code{TRUE}}
#'
#' @param debug \code{TRUE} if code debugging messages must be printed.
#' \cr \emph{Default: \code{FALSE}}
#'
#' @details At each regularization parater value \code{lambda}, cyclic coordinate
#' descent is used to update the parameters until convergence. The intercept and
#' the scale parameter are never regularized.
#' The obtained solution is used to initialize the parameters at the next \code{
#' lambda} value.
#' @return Returns a S3 object iregnet with the following elements:\cr
#' \tabular{ll}{
#'  \code{beta} \tab Matrix of size \code{(n_vars+1) * num_lambda} containing
#'  intercept, coefficients of \code{X} for each \code{lambda} in the fit model.
#'    \cr
#'  \code{call} \tab Copy of the call that produced this object. \cr
#'  \code{lambda} \tab Vector of size \code{num_lambda} of (calculated or
#'   supplied) regularization parameter \code{lambda} values. \cr
#'  \code{loglik} \tab Vector of size \code{num_lambda} of log-likelihoods of
#'    the fit at each \code{lambda} value, excluding the contribution of the
#'    penalty terms. \cr
#'  \code{num_lambda} \tab Number of \code{lambda} values. \cr
#'  \code{n_iters} \tab Vector of size \code{num_lambda} of number of iterations
#'    taken at each \code{lambda}. \cr
#'  \code{scale} \tab Vector of size \code{num_lambda} of estimated
#'    scale at each \code{lambda} value, if \code{estimate_scale == TRUE}. Same as
#'    \code{scale_init} otherwise. \cr
#'  \code{scale_init} \tab Initial value (calculated or supplied) of \code{scale}. \cr
#'  \code{estimate_scale} \tab \code{TRUE} if the \code{scale} was estimated. \cr
#'  \code{error_status} \tab The error status. \code{0} denotes no errors.
#'    \code{-1} denotes that convergence was not reached in \code{maxiter}. \cr
#' }
#' @author
#' Anuj Khare, Toby Dylan Hocking, Jelle Goeman. \cr
#' Maintainer: Anuj Khare \email{khareanuj18@gmail.com}
#'
#' @section References:
#' Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#' Generalized Linear Models via Coordinate Descent,
#' \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}
#'
#' Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization
#' Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal
#' of Statistical Software, Vol. 39(5) 1-13
#' \url{http://www.jstatsoft.org/v39/i05/}
#'
#' @seealso
#' \code{\link{predict.iregnet}}, \code{cv.iregnet}, \code{\link{plot.iregnet}}
#'
#' @examples
#' # y can be a 2 column matrix.
#' set.seed(10)
#' X <- matrix(rnorm(50), 10, 5)
#' y <- matrix(rnorm(20), 10, 2)
#' y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
#' fit1 <- iregnet(X, y)
#'
#' # Surv objects from survival are also supported.
#' data("ovarian")
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit2 <- iregnet(X, y)
#'
#' # Log-Gaussian is same as Gaussian with log-transformed data
#' set.seed(10)
#' X <- matrix(rnorm(50), 10, 5)
#' y <- matrix(abs(rnorm(20)), 10, 2)
#' y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
#' fit3 <- iregnet(X, log(y), "gaussian")
#' fit4 <- iregnet(X, y, "loggaussian")
#'
#' # Scale parameter can be fixed by setting the estimate_scale flag.
#' set.seed(10)
#' X <- matrix(rnorm(50), 10, 5)
#' y <- matrix(rnorm(20), 10, 2)
#' y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
#' fit5 <- iregnet(X, y, scale_init=1, estimate_scale=FALSE)
#'
iregnet <- function(x, y,
                    family=c("gaussian", "logistic", "loggaussian", "loglogistic", "extreme_value", "exponential", "weibull"),
                    alpha=1, lambda=double(0), num_lambda=100, intercept=TRUE, standardize=TRUE, scale_init=NA, estimate_scale=TRUE,
                    maxiter=1e3, threshold=1e-4, unreg_sol=TRUE, eps_lambda=NA, debug=0) {

  # Parameter validation ===============================================
  stopifnot_error("alpha should be between 0 and 1", 0 <= alpha, alpha <= 1)
  stopifnot_error("num_lambda > 0 is required", num_lambda > 0)
  stopifnot_error("lambdas must be numeric", is.numeric(lambda))
  stopifnot_error("intercept must be a boolean flag", is.logical(intercept))
  stopifnot_error("standardize must be a boolean flag", is.logical(standardize))
  stopifnot_error("unreg_sol must be a boolean flag", is.logical(unreg_sol))
  stopifnot_error("estimate_scale must be a boolean flag", is.logical(estimate_scale))
  stopifnot_error("threshold must be positive", threshold > 0)

  # if (estimate_scale == FALSE && is.na(scale_init))
  #   stop("Value of scale required if scale is not estimated")

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
    stopifnot_error("y should be a 2 column matrix", ncol(y) == 2)
  }
  stopifnot_error("nrow(y) = nrow(x) is not true", nrow(y) == n_obs)

  temp <- y[0]; y[0] <- 1; y[0] <- temp # FIXME: We need deep copy of y, otherwise C++ modifies it
  stopifnot_error("y should be positive for the given family",
                  !(family %in% names(transformed_distributions) && any(y[!is.na(y)]<0)))

  # Fix scale for exponential: (least) extreme value distribution with scale = 1
  if (family == "exponential") {
    cat("Exponential distribution: fixing scale to 1\n")
    estimate_scale <- FALSE
    scale_init <- 1
  }
  # Transform the outputs, and get new dist
  if (family %in% names(transformed_distributions)) {
    trans <- transformed_distributions[[family]]
    y <- trans$trans(y)
    family <- trans$dist
  }

  # Get column names
  varnames <- colnames(x)
  if (is.null(varnames)) {
    varnames <- paste('x', 1: n_vars, sep='')
  }

  # Append col of 1's for the intercept
  if (intercept) {
    x <- cbind(rep(1, n_obs), x)
    varnames = c("(Intercept)", varnames)
  }

  if (is.na(eps_lambda))
    eps_lambda <- ifelse(n_obs < n_vars, 0.01, 0.0001)
  stopifnot_error("eps_lambda should be between 0 and 1", 0 <= eps_lambda && eps_lambda < 1)

  # Call the actual fit method
  fit <- fit_cpp(x, y, family, alpha, lambda_path=lambda, num_lambda=num_lambda, intercept=intercept,
                 out_status=status, scale_init=scale_init, max_iter=maxiter, threshold=threshold,
                 flag_standardize_x=standardize, estimate_scale=estimate_scale, unreg_sol=unreg_sol,
                 eps_lambda=eps_lambda, debug=debug);

  fit$call <- match.call()
  fit$intercept <- intercept
  fit$family <- family
  rownames(fit$beta) <- varnames
  class(fit) <- 'iregnet'
  fit
}
