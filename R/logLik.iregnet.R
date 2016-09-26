#' Extract Log-Likelihood
#' 
#' @method
logLik.iregnet <- function(object, ...) {
  stopifnot_error("Invalid / no object provided",
                  class(object) == "iregnet")

  object$loglik
}
