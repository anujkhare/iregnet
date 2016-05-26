library("testthat")

# Tests for the fit_cpp C++ function
# Note: the debug flag values are given by the enum in iregnet.h

y <- cbind(c(NA, 1, 2, NA, 4), c(NA, NA, 3, 3, 4));
x <- cbind(c(NA, 1, 2, NA, 4), c(1, NA, 3, 3, 4));
l = fit_cpp(x, y, "gaussian", 1, flag_debug = 3);

test_that("Censoring types are calculated correctly", {
  expect_equal(l$error_status, 0);
  expect_equal(l$censoring_types, c(4, 0, 3, 2, 1));
})
