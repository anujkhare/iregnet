#' Extract model coefficients from iregnet object
#'
#' @description
#' Returns the model coefficients for the fit \link{iregnet} object. One model
#' is fit for each lambda value in \code{object$lambda}. Each of these is
#' returned in a single matrix of size \code{nvars} * \code{nobs}.
#'
#' @param object The S3 object of type \code{iregnet} returned by the
#' \code{iregnet} method.
#' @param ... Optional arguments. Currently unused.
#'
#' @method coef iregnet
coef.iregnet <- function(object, ...) {
  stopifnot_error("Invalid / no object provided",
                  class(object) == "iregnet")

  object$beta
}
