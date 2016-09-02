# Print the iregnet fit
# <TODO>
print.iregnet <- function(x) {
  cat('\nCall:', deparse(x$call), '\n')
  out_matrix <- t(with(x, rbind(niters=n_iters, lambda=lambda, scale=scale, beta)))
  print(out_matrix[1:10, ])
}
