# iregnet funciton
# USAGE
# Authors
# Borrowed parts from survival survreg and glmnet

# TODO: optimization params need to be implemented
# TODO: checks like - can't provide scale with exp
iregnet <- function(x, y,
                    family=c("gaussian", "logistic", "loggaussian", "extreme value", "exponential"),
                    alpha=1, flag_debug=0, intercept=T, standardize=F, scale=NA) {

  # Parameter validation ===============================================
  # alpha should be between 0 and 1
  if (alpha > 1) {
    warning ("alpha > 1; setting to 1");
    alpha <- 1;
  }
  if (alpha < 0) {
    warning ("alpha<0; set to 0");
    alpha <- 0;
  }
  alpha <- as.double(alpha);

  # family should be one of ""
  family <- match.arg(family)
  # print(family)

  # x should be a matrix with atleast two columns (TODO: why glmnet requires 2 or more cols?)
    # TODO: NAs in x? should be numeric
  dim_x <- dim(x);    # LIST
  if (is.null(dim_x) | dim_x[2] <= 1)
    stop("x should be a matrix with 2 or more columns");

  n_obs  <- dim_x[1];   # LIST
  n_vars <- dim_x[2];   # LIST

  # y should be a matrix with 2 columns correspoding to the left and right times
  # NA or Inf/-Inf can be used for censored entries
  dim_y <- dim(y);      # LIST
  if (is.null(dim_y) | dim_y[2] != 2)
   stop("y should be a matrix with 2 columns");

  # length of x and y should be the same
  if (dim_y[1] != n_obs)
    stop("number of rows in x and y don't match");

  # Call the actual fit method
  #fit_cpp(x, y, family, alpha, intercept, standardize, scale=scale);
  # fit_cpp(cbind(rep(1, length(y)), x), y, family, alpha, intercept, standardize, scale=scale);
  if (intercept) {
    fit_cpp(cbind(rep(1, n_obs), x), y, family, alpha, intercept=TRUE, scale=scale,
            standardize = standardize);
  } else {
    fit_cpp(x, y, family, alpha, intercept=FALSE, scale=scale,
            standardize = standardize);
  }
}
