library(iregnet)

test_that("Output y is validated properly", {
  data("ovarian", package="survival")
  x <- cbind(ovarian$ecog.ps, ovarian$rx)
  y_l <- ovarian$futime
  y_r <- ovarian$futime
  y_r[ovarian$fustat == 0] <- NA

  expect_error(iregnet(x, 10), "y should be a 2 column matrix, or a Surv object")
  # 2 column matrix
  expect_error(iregnet(x, cbind(y_l, y_r-1)), "Invalid interval*")
  expect_error(iregnet(x, cbind(y_r, y_r)), "Invalid interval*")
  expect_error(iregnet(x, cbind(y_l, y_l, y_l)), "y should be a 2 column matrix")

  # Surv object
  expect_error(iregnet(x, Surv(y_l, ovarian$fustat, type="mstate")), "Unsupported censoring type from Surv")

  # Size wrt X
  expect_error(iregnet(x, cbind(y_l[1:2], y_r[1:2])), "*")
  expect_error(iregnet(x, Surv(y_l[1:2], y_r[1:2], type="interval2")), "*")

  # postivity
  expect_error(iregnet(x, cbind(-y_l, y_r), family="loglogistic"), "y should be positive for the given family")
});

test_that("Fit parameters are validated properly", {
  data("ovarian", package="survival")
  x <- cbind(ovarian$ecog.ps, ovarian$rx)
  y <- Surv(ovarian$futime, ovarian$fustat)

  # lambda
  expect_error(iregnet(x, y, lambda='s'), "lambdas must be numeric")
  expect_error(iregnet(x, y, lambda=c(0.1, -0.3)), "lambdas must be positive and decreasing")
  expect_error(iregnet(x, y, lambda=c(0.1, 0.3)), "lambdas must be positive and decreasing")

  # misc
  expect_error(iregnet(x, y, alpha=-1), "alpha should be between 0 and 1")
  expect_error(iregnet(x, y, eps_lambda=-1), "eps_lambda should be between 0 and 1")
});
