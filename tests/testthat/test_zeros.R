library(iregnet)

source('get_xy.R')

test_zeros <- function(x, y, dist="gaussian") {
  fit_i <- iregnet(x, y, family=dist,
                   alpha=1, intercept = T, threshold=1e-4, debug=0)

  # TODO: Currently, the first solution is for the initial fit. remove it!
  first_solution  <- as.double(fit_i$beta[, 2])
  second_solution <- as.double(fit_i$beta[, 3])
  nvars <- length(first_solution)

  expect_equal(first_solution[2: nvars], rep(0, nvars - 1))
  # want at least one non-zero coefficient besides intercept
  expect_equal(any(second_solution[2: nvars] != 0), T)
}

test_that("G right censored: First solution is intercept only, rest are not", {
  data("ovarian", package="survival")
  x <- cbind(ovarian$ecog.ps, ovarian$rx)
  y <- Surv(ovarian$futime, ovarian$fustat)
  test_zeros(x, y)
})

test_that("G exact: first solution is intercept only, rest are not", {
  set.seed(10)
  xy <- get_xy(30, 6, "none", standardize=F)
  test_zeros(xy$x, xy$surv)
})

test_that("G left: first solution is intercept only, rest are not", {
  set.seed(10)
  for (nvars in 6:10) {
    xy <- get_xy(30, nvars, "left", standardize=F)
    test_zeros(xy$x, xy$surv)
  }
})
test_that("G right: first solution is intercept only, rest are not", {
  set.seed(10)
  xy <- get_xy(30, 6, "right", standardize=F)
  test_zeros(xy$x, xy$surv)
})

test_that("exact: first solution is intercept only, rest are not", {
  set.seed(10)
  xy <- get_xy(30, 6, "none", standardize=F, positive=T)
  test_zeros(xy$x, xy$surv, "logistic")
  test_zeros(xy$x, xy$surv, "loglogistic")
  test_zeros(xy$x, xy$surv, "exponential")
})
