library("iregnet")
library("survival")
library("glmnet")

# Right censored data from Survival, get into our format first
data("ovarian")
y_l <- y_r <- ovarian$futime
y_r[ovarian$fustat == 0] <- NA
y <- cbind(y_l, y_r)

x <- cbind(ovarian$ecog.ps, ovarian$rx)

fit_s <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, dist = "gaussian")
fit_i <- iregnet(x, y, family="gaussian", alpha=1, intercept = T, scale = fit_s$scale)
# print (as.double(fit_s$coefficients))
# print (fit_i$beta[, fit_i$num_lambda + 1])
test_that("iregnet calculates correct coefficients", {
  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
})

# TODO: integrate Surv support into iregnet, for now generate eq. dists:
get_xy <- function(n_obs, n_vars, type = "right") {
  x <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars)
  y <- rnorm(n_obs)

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

test_that("Gaussian, left censored data - coefficients are calculated correctly:", {
  set.seed(55)

  # n_obs >> n_vars
  for (n_vars in 2:10) {
    xy <- get_xy(40, n_vars, "left")
    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
    # print(fit_s$coefficients)
  }

  # n_obs >= n_vars, but smaller - TODO: FAILING!
  xy <- get_xy(11, n_vars, "left")
  fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian", control = survreg.control(iter.max=1000))
  fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
})

test_that("Gaussian, right censored data - coefficients are calculated correctly:", {
  set.seed(55)

  # n_obs >> n_vars
  for (n_vars in 2:10) {
    xy <- get_xy(40, n_vars, "right")
    fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
    fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
    expect_equal(as.double(fit_s$coefficients),
                 fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
    # print(fit_s$coefficients)
  }

  # n_obs >= n_vars, but smaller - TODO: FAILING!
  xy <- get_xy(11, n_vars, "right")
  fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian", control = survreg.control(maxiter=1000, iter.max=1000))
  fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
})


test_that("Gaussian, exact data - coefficients are calculated correctly:", {
  set.seed(115)

  n_vars <- 5;
  xy <- get_xy(30, n_vars, "none")

  fit_s <- survreg(xy$surv ~ xy$x, dist = "gaussian")
  fit_i <- iregnet(xy$x, xy$y, alpha = 1, intercept = T, scale = fit_s$scale)
  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)
  expect_equal(as.double(fit_s$coefficients),
               fit_i$beta[, fit_i$num_lambda + 1], tolerance = 1e-3)

  # TODO: VERY different from glmnet
  fit_g <- glmnet(xy$x, xy$y[, 1], "gaussian")
  fit_g <- glmnet(xy$x, xy$y[, 1], "gaussian")
  print (fit_i$beta)
  print (coef(fit_g))
})
