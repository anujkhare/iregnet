#' @title Print the iregnet fit
#' @export
#' @description
#' Prints a summary of the results of the 
#' @import utils
#' @param x The result of an iregnet fit.
#'
#' @param ... Optional parameters for print. If \code{n} is supplied, only the
#' first \code{n} models will be printed.
#' @method print iregnet
print.iregnet <- function(x, ...) {
  stopifnot_error("Invalid / no x object provided", !missing(x),
                  class(x) == "iregnet")
  cat('\nCall:', deparse(x$call), '\n\n')
  varargs <- list(...)
  n <- varargs$n
  if(is.null(n))
    n <- x$num_lambda
  print(head(with(x, cbind(lambda, scale, loglik, n_iters, t(beta))), ...))
}
