library("testthat")

# Tests for the fit_cpp C++ function
# Note: the debug flag values are given by the enum in iregnet.h

y <- cbind(c(NA, 1, 2, NA, 4), c(NA, NA, 3, 3, 4));
x <- cbind(c(NA, 1, 2, NA, 4), c(1, NA, 3, 3, 4));
l = fit_cpp(x, y, "gaussian", 1, flag_debug = 3);

# TODO: validations for y
#test_that("Censoring types are calculated correctly", {
#  expect_equal(l$error_status, 0);
#  expect_equal(l$censoring_types, c(4, 0, 3, 2, 1));
#})

# Test internal variables using ovarian data
data("ovarian")
# test_that("Right censored data... ")
# y_l <- y_r <- ovarian$futime
# y_r[ovarian$fustat == 0] <- NA
# y <- cbind(y_l, y_r)
# X <- cbind(ovarian$ecog.ps, ovarian$rx)

# iregnet(X, y, family="exponential")

# test_that("No censoring waala data", {})
# X <- cbind(ovarian$ecog.ps, ovarian$rx)
# y <- ovarian$futime / sqrt(var(ovarian$futime))
# y <- cbind(y, y)
# print (y)

load("/home/kcm/code/xy")
X <- a[[1]]
y <- a[[2]]
#y <- y / as.double(sqrt(var(y)))
y <- cbind(y, y)
fit <- iregnet(X, y, family = "gaussian", alpha=1, scale = 1, intercept = T)
print (fit$lambda)
print (fit$n_iters)
print (fit$beta)

test_that("Derivatives wrt eta are calculated correctly", {

  y_l <- c(4.07754, 4.74493, 5.04986, 6.04263, 6.06611, 6.10479, 6.13988, 6.16331, 6.16752, 6.33328, 6.45834, 6.61204, 6.64509, 6.64639, 6.68835, 6.7511, 6.94698, 7.00851, 7.02909, 7.09506, 7.11233, 5.59099, 5.79606, 5.86647, 5.8999, 5.93225)
  y_r <- c(4.07754, 4.74493, 5.04986, NA, 6.06611, NA, 6.13988, 6.16331, NA, 6.33328, 6.45834, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 5.59099, 5.79606, 5.86647, 5.8999, NA)

  censoring <- c(1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0)

  # Values taken from some random iteration of Survival survreg
  eta <- c(6.462, 6.462, 6.370, 7.082, 6.462, 6.370, 6.990, 6.990, 6.462, 6.990, 6.370, 7.082, 6.990, 7.082, 6.462, 6.370, 6.370, 6.462, 7.082, 7.082, 6.990, 6.370, 6.462, 6.990, 7.082, 7.082)
  mu <-  c(-0.908, -0.820, -0.733, 0.354, -0.327, 0.767, -0.573, -0.563, 0.745, -0.482, 0.092,  0.625,  0.708,  0.647,  1.254,  1.464,  1.781,  1.728,  0.948,  1.013,  1.130,  -0.541, -0.486, -0.675, -0.693, 0.317)
  w <- c(-0.092, -0.180, -0.267, -0.354, -0.673, -0.767, -0.427, -0.437, -0.745, -0.518, -1.092, -0.625, -0.708, -0.647, -1.254, -1.464, -1.781, -1.728, -0.948, -1.013, -1.130, -0.459, -0.514, -0.325, -0.307, -0.317)
  fit <- compute_grad_response_cpp(y_l, y_r, eta, 1, censoring, "extreme_value")

  expect_equal(fit$mu, mu, tolerance=1e-3)   # I copied values rounded upto 3 decimals
  expect_equal(fit$w, w, tolerance=1e-3)

  #print(cbind(w, fit$z));
})
#print (cbind(fit$mu, mu))
#print (cbind(fit$w, w_1))
#print (cbind(fit$z, z))

test_that("Densities are calculated correctly", {

  # Gaussian
  x <- -10:10

  df <- exp(-x * x / 2) * (-x) / sqrt(2*pi)               # derivatives
  calculated <- compute_densities(x, 2, "gaussian")       # F, 1-F, f, f'
  expected <- cbind(pnorm(x), 1-pnorm(x), dnorm(x), df)
  # expect_equal(compute_densities(x, 2, "gaussian"),
  #              cbind(pnorm(x), 1-pnorm(x), dnorm(x), df)
  #              )
  expect_equal(calculated[, 1], expected[, 1])
  expect_equal(calculated[, 2], expected[, 2])
  expect_equal(calculated[, 3], expected[, 3])
  expect_equal(calculated[, 4], expected[, 4])

  # # Lognormal
  # x <- 1:10
  # all.equal(dsurvreg(x, 1, 5, 'lognormal'), dlnorm(x, 1, 5))
  # all.equal(psurvreg(x, 1, 5, 'lognormal'), plnorm(x, 1, 5))
  # all.equal(qsurvreg(p, 1, 5, 'lognormal'), qlnorm(p, 1, 5))
  #
})
# # Weibull
# lambda <- exp(-2)
# rho    <- 1/3
# temp <- (lambda*x)^rho
# all.equal(psurvreg(x, 2, 3), 1- exp(-temp))
# all.equal(dsurvreg(x, 2, 3), lambda*rho*(lambda*x)^(rho-1)*exp(-temp))
#
