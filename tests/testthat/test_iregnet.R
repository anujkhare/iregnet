library("iregnet")
library("survival")
library("glmnet")

# TODO: You should check if the suggested packages are present, and only then use them for testing

std <- F
get_xy <- function(n_obs, n_vars, type = "right", standardize=std) {
  x <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars)
  y <- rnorm(n_obs)

  # standardize x and y
  if (standardize == T) {
    for (i in 1:ncol(x)) {
      x[, i] <- (x[, i] - mean(x[, i])) / sd(x[, i]);
    }

    y <- (y - mean(y))
    y <- y / sd(y)
  }

  if (type == "none") {
    status = rep(1, length(y))
    y_surv <- Surv(time = y, event = status, type = "right")

  } else {
    status <- sample(c(0, 1), size=n_obs, replace=T)
    y_surv <- Surv(time = y, event = status, type = type)
  }

  # get the y matrix
  y <- cbind(y, y)
  if (type=="right") {
    y[status == 0, 2] = NA
  } else {
    y[status == 0, 1] = NA
  }

  return (list("x" = x, "y" = y, "surv" = y_surv))
}

test_that("iregnet calculates correct coefficients for ovarian data wrt survival", {
  data("ovarian")
  x <- cbind(ovarian$ecog.ps, ovarian$rx)

  fit_s <- survreg(Surv(futime, fustat) ~ x, data = ovarian, dist = "gaussian")
  fit_i <- iregnet(x, Surv(ovarian$futime, ovarian$fustat), family="gaussian", alpha=1, intercept = T, threshold=1e-4)

  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
})


test_that("Gaussian, left censored data - coefficients are calculated correctly wrt survival:", {
  set.seed(55)

  # n_obs >> n_vars
  for (n_vars in 2:10) {
  # n_vars <- 5
    xy <- get_xy(40, n_vars, "left", standardize=std)
    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$surv, "gaussian", alpha = 1, intercept = T)
    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
  }

  # n_obs >= n_vars, but smaller - TODO: FAILING!
  # xy <- get_xy(11, n_vars, "left")
  # fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian", control = survreg.control(iter.max=1000))
  # fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
  # expect_equal(as.double(fit_s$coefficients),
  #              fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
})

test_that("Gaussian, right censored data - coefficients are calculated correctly wrt survival:", {
  set.seed(55)

  # n_obs >> n_vars
  for (n_vars in 2:10) {
    xy <- get_xy(40, n_vars, "right", standardize=std)
    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    # fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
    fit_i <- iregnet(xy$x, xy$y, "gaussian", alpha = 1, intercept = T)
    # print(fit_i)
    # print(fit_s)
    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
  }

  # n_obs >= n_vars, but smaller - TODO: FAILING!
  # xy <- get_xy(11, n_vars, "right")
  # fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian", control = survreg.control(maxiter=1000, iter.max=1000))
  # fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
  # expect_equal(as.double(fit_s$coefficients),
  #              fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
})


test_that("ElemStatsLearn data - coefficients are calculated correctly wrt survival and glmnet:", {
  alpha <- 0.6

  data(prostate,package="ElemStatLearn")
  pros <- subset(prostate,select=-train,train==TRUE)
  ycol <- which(names(pros)=="lpsa")
  X.unscaled <- as.matrix(pros[-ycol])
  y.unscaled <- pros[[ycol]]
  M <- matrix(
    colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
  X.centered <- X.unscaled - M
  sd.vec <- apply(X.unscaled, 2, sd)
  S <- diag(1/sd.vec)
  X.scaled <- X.centered %*% S
  dimnames(X.scaled) <- dimnames(X.unscaled)
  m <- mean(y.unscaled)
  sigma <- sd(y.unscaled)
  y.scaled <- (y.unscaled - m)/sigma

  X <- X.scaled
  y <- y.scaled

  fit_s <- survreg(Surv(y, rep(1, length(y))) ~ X, dist = "gaussian")
  fit_i <- iregnet(X, cbind(y, y), "gaussian", maxiter=1e5, thresh=1e-7, standardize=F, alpha=alpha, scale=1, estimate_scale=F)

  lambda_path <- fit_i$lambda * (fit_i$scale ** 2)
  # lambda_path <- fit_i$lambda

  fit_g <- glmnet(X, y, "gaussian", lambda = lambda_path, standardize=F, maxit=1e5, thresh=1e-7, alpha=alpha)

  expect_equal(as.double(fit_s$coefficients), fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
  expect_equal(as.double(fit_i$beta), as.double(coef(fit_g)), tolerance=1e-3)
})


test_that("Gaussian, exact data - coefficients are calculated correctly wrt survival and glmnet:", {
  set.seed(115)

  n_vars <- 5
  # for (n_vars in 5:10)
  {
    # FIXME: to get same results from Glmnet and Iregnet, y should have 0 mean and 1 var
    xy <- get_xy(30, n_vars, "none", standardize=std)

    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$y, "gaussian", alpha = 1, intercept = T, thresh=1e-4, standardize=T)

    lambda_path <- fit_i$lambda * (fit_i$scale ** 2)
    fit_g <- glmnet(xy$x, xy$y[, 1], "gaussian", lambda=lambda_path, thresh=1e-7, standardize=T)

    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)

    expect_equal(as.double(fit_i$beta), as.double(coef(fit_g)), tolerance=1e-3)
  }
})

test_that("Extreme value dist - coefficients are calculated correctly wrt survival:", {
  set.seed(115)

  n_vars <- 5
  # for (n_vars in 5:10)
  {
    # FIXME: interval censored not tested
    xy <- get_xy(30, n_vars, "right", standardize=std)

    fit_s <- survreg(xy$surv ~ xy$x, dist = "extreme")
    fit_i <- iregnet(xy$x, xy$y, "extreme", alpha = 1, intercept = T, thresh=1e-4, standardize=T)

    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
  }

    xy <- get_xy(30, n_vars, "left", standardize=std)

    fit_s <- survreg(xy$surv ~ xy$x, dist = "extreme")
    fit_i <- iregnet(xy$x, xy$y, "extreme", alpha = 1, intercept = T, thresh=1e-4, standardize=T)

    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda], tolerance = 1e-3)
})
