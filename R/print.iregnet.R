# Print the iregnet fit
# <TODO>
#' @method
print.iregnet <- function(x, ...) {
  stopifnot_error("Invalid / no x object provided", !missing(x),
                  class(x) == "iregnet")
  cat('\nCall:', deparse(x$call), '\n\n')
  tidy_matrix <- tidymatrix(x)
  print(head(tidy_matrix, n=10))
}
