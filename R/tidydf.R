tidydf <- function(fit) {
  stopifnot_error("Invalid / no fit object provided",
                  class(fit) == "iregnet")

  # Don't include intercept in norm (arclength) since it is never regularized.
  start_index <- as.integer(fit$intercept) + 1
  n <- nrow(fit$beta)
  arclength <- apply(fit$beta, 2, function(x) sum(abs(x[start_index: n])))
  tidy.df <- with(fit, data.frame(
    weight=as.numeric(t(beta)),
    lambda,
    arclength,
    variable=rownames(beta)[as.integer(col(t(beta)))]
    ))
  tidy.df
}
