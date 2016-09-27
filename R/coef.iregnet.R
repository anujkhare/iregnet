#' @method coef iregnet
coef.iregnet <- function(object, ...) {
  stopifnot_error("Invalid / no object provided",
                  class(object) == "iregnet")

  object$beta
}
