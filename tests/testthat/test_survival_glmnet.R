library(iregnet)
library(survival)
library(glmnet)

# TODO: You should check if the suggested packages are present, and only then use them for testing
std <- F
get_xy <- function(n_obs, n_vars, type = c('right', 'left', 'none', 'interval'), standardize=std, positive=F) {
  type <- match.arg(type)

  x <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars)
  y <- rnorm(n_obs)
  y_u <- rnorm(n_obs)

  # standardize x and y
  if (standardize == T) {
    for (i in 1:ncol(x)) {
      x[, i] <- (x[, i] - mean(x[, i])) / sd(x[, i]);
    }

    y <- (y - mean(y))
    y <- y / sd(y)

    if (type == 'interval') {
      y_u <- (y_u - mean(y_u))
      y_u <-  y_u / sd(y_u)
    }
  }

  if (positive) {
    y <- abs(y)
    y_u <- abs(y_u)
  }

  if (type == "none") {
    status = rep(1, length(y))
    y_surv <- Surv(time = y, event = status, type = "right")

  } else if (type == 'interval') {
    status <- sample(c(0, 1, 2, 3), size=n_obs, replace=T)
    y <- cbind(y, y_u)
    y <- t(apply(y, 1, sort))    # make sure intervals are increasing in time
    y_surv <- Surv(time=y[, 1], time2=y[, 2], event = status, type = 'interval')

  } else {    # left or right censored
    status <- sample(c(0, 1), size=n_obs, replace=T)
    y_surv <- Surv(time = y, event = status, type = type)
  }

  # get the y matrix
  y <- cbind(y, y)
  if (type=="right") {
    y[status == 0, 2] = NA
  } else if (type == 'interval') {
    y <- NA  # NOTE: Not implemented, not needed, use the Surv object!
  } else {
    y[status == 0, 1] = NA
  }

  return (list("x" = x, "y" = y, "surv" = y_surv))
}

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

  #for (n_vars in 5:10)
  n_vars <- 5
  {
    xy <- get_xy(30, n_vars, "none", standardize=std)
    xy1 <- xy

    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$y, "gaussian", alpha = 1, intercept = T, thresh=1e-4, standardize=T,
    debug=F)

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
