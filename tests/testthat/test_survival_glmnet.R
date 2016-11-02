context("\nComparision with survival and glmnet")
library(iregnet)
library(survival)
library(glmnet)

source('get_xy.R')

# TODO: You should check if the suggested packages are present, and only then use them for testing
std <- F

test_that("survival::ovarian data: iregnet calculates correct coefficients wrt survival", {
  data("ovarian")
  x <- cbind(ovarian$ecog.ps, ovarian$rx)

  fit_s <- survreg(Surv(futime, fustat) ~ x, data = ovarian, dist = "gaussian")
  fit_i <- iregnet(x, Surv(ovarian$futime, ovarian$fustat), family="gaussian", alpha=1, intercept = T, threshold=1e-4)

  expect_equal(fit_s$coefficients,
               fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
})

test_that("Gaussian, exact data - coefficients are calculated correctly wrt survival and glmnet:", {
  set.seed(115)

  for (n_vars in 5:10)
  {
    xy <- get_xy(30, n_vars, "none", standardize=std)

    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$y, "gaussian", alpha = 1, intercept = T, thresh=1e-4, standardize=T,
    debug=F)

    # Note: fit$lambda returned by iregnet are scaled by (1/scale^2) in the
    # case of gaussian exact data. Rescale to compare the solutions.
    lambda_path <- fit_i$lambda * (fit_i$scale ** 2)
    fit_g <- glmnet(xy$x, xy$y[, 1], "gaussian", lambda=lambda_path, thresh=1e-7, standardize=T)

    expect_equal(as.double(fit_s$coefficients),
                 as.double(fit_i$beta[, fit_i$num_lambda]), tolerance = 1e-3)

    expect_equal(as.double(fit_i$beta), as.double(coef(fit_g)), tolerance=1e-3)
  }
})

test_wrt_survival <- function(dist, types = c('left', 'right', 'none'), n_vars_list = 5:10,
                              std=F, nobs = 30, seed_test = 185, tol = 1e-3)
{
  set.seed(seed_test)
  for (i in seq_along(n_vars_list))
  {
    for (type in types) {
      xy <- get_xy(nobs, n_vars_list[i], type, standardize=std, T)

      fit_s <- survreg(xy$surv ~ xy$x, dist = dist)
      fit_i <- iregnet(xy$x, xy$surv, dist, alpha = 1, intercept = T, thresh=1e-7, standardize=T)

      expect_equal(as.double(fit_s$coefficients),
                   as.double(fit_i$beta[, fit_i$num_lambda]), tolerance = tol)
    }
  }
}

test_that("Gaussian dist - unregularized coefficients are calculated correctly wrt survival:", {
   test_wrt_survival("gaussian", c('left', 'right', 'interval'), 2:10)
})

test_that("Logistic dist - unregularized coefficients are calculated correctly wrt survival:", {
   test_wrt_survival("logistic", c('left', 'right', 'interval'), 2:10)
})

test_that("LogGaussian - coefficients are calculated correctly wrt survival:", {
   test_wrt_survival("loggaussian", c('left', 'right', 'interval'), 2:10)
})

test_that("LogLogistic - coefficients are calculated correctly wrt survival:", {
  test_wrt_survival("loglogistic", c('left', 'right', 'interval'), 2:10)
})

test_that("Exponential dist - unregularized coefficients are calculated correctly wrt survival:", {
   test_wrt_survival("exponential", c('left', 'right', 'interval'), 2:10)
})
