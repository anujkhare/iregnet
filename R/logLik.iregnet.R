#' @title Extract Log-Likelihood
#'
#' @description
#' Extracts the Log-Likelihood from an object of class "iregnet" which has
#' already been fit.
#' 
#' @param object The result of an iregnet fit.
#'
#' @param ... Optional parameters, \emph{unused}.
#' @method logLik iregnet
logLik.iregnet <- function(object, ...) {
  stopifnot_error("Invalid / no object provided",
                  class(object) == "iregnet")

  object$loglik
}
