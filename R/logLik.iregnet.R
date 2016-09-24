logLik.iregnet <- function(fit, ...) {
  stopifnot_error("Invalid / no fit object provided",
                  class(fit) == "iregnet")

  fit$loglik
}
