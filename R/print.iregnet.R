# Print the iregnet fit
# <TODO>
print.iregnet <- function(fit) {
  stopifnot_error("Invalid / no fit object provided", !missing(fit),
                  class(fit) == "iregnet")
  cat('\nCall:', deparse(fit$call), '\n\n')
  tidy_matrix <- tidymatrix(fit)
  print(head(tidy_matrix, n=10))
}
