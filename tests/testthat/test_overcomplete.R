library(iregnet)
library(survival)

source('get_xy.R')

test_that("Gaussian exact - nvars > nobs", {
  xy <- get_xy(10, 15, "none", standardize=F)

  fit_i <- iregnet(xy$x, xy$y, "gaussian", alpha = 1, intercept = T,
                   thresh=1e-4, standardize=T, debug=F)

  lambda_path <- fit_i$lambda * (fit_i$scale ** 2)
  fit_g <- glmnet(xy$x, xy$y[, 1], "gaussian", lambda=lambda_path,
                  thresh=1e-7, standardize=T)

  expect_equal(as.double(fit_i$beta), as.double(coef(fit_g)), tolerance=1e-3)
})
