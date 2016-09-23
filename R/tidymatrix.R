tidymatrix <- function(fit) {
  stopifnot_error("Invalid / no fit object provided",
                  class(fit) == "iregnet")

  # Don't include intercept in norm (arclength) since it is never regularized.
  start_index <- as.integer(fit$intercept) + 1
  n <- nrow(fit$beta)
  norm <- apply(fit$beta, 2, function(x) sum(abs(x[start_index: n])))
  tidy_matrix <- t(with(fit, rbind(lambda, n_iters, scale, loglik, norm, coef=beta)))
  tidy_matrix
}
