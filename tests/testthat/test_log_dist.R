context("\nLog-dists")
library(iregnet)

test_that("Log* with y is same as * with log(y)", {
  data(ovarian)
  X <- cbind(ovarian$ecog.ps, ovarian$rx)
  y <- Surv(ovarian$futime, ovarian$fustat)
  y_log <- Surv(log(ovarian$futime), ovarian$fustat)
  thresh <- 1e-7

  dists <- c("gaussian", "logistic")
  log_dists <- c("loggaussian", "loglogistic")
  for (i in seq_along(dists)) {
    fit_il <- iregnet(X, y, log_dists[i], thresh=thresh)
    fit_i <- iregnet(X, y_log, dists[i], thresh=thresh)
    fit_s <- survreg(y_log ~ X, dist=dists[i])
    expect_equal(fit_il$beta, fit_i$beta, tolerance=1e-3)
    expect_equal(as.double(fit_s$coefficients),
                 as.double(fit_i$beta[, fit_i$num_lambda]), tolerance = 1e-3)
  }
})
