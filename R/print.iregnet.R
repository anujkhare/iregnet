#' @title Print the iregnet fit
#'
#' @description
#' Prints a summary of the results of the 
#'
#' @param object The result of an iregnet fit.
#'
#' @param ... Optional parameters for print.
#' @method print iregnet
print.iregnet <- function(x, ...) {
  stopifnot_error("Invalid / no x object provided", !missing(x),
                  class(x) == "iregnet")
  cat('\nCall:', deparse(x$call), '\n\n')
  tidy_matrix <- tidymatrix(x)
  print(tidy_matrix, ...)
}
