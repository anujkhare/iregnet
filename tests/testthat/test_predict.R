context("\nPredict function")
library("iregnet")
library("testthat")
test_that("predict.iregnet returns expected results", {
  x <- matrix(rnorm(200), 20, 10)
  y <- rnorm(20)

  fit_i <- iregnet(x, cbind(y, y), "gaussian", standardize=T)

  newx <- matrix(rnorm(80), 8, 10)
  expect_error(predict(fit_i, newx, 1e11), "Lambda values must be those used in fit")
  expect_error(predict(fit_i, cbind(newx, newx)),"features missing but needed for prediction: ")

  #inds <- c(2, 10, 30, 40)
  #pred_i <- predict(fit_i, newx, fit_i$lambda[inds], type="response")
  #newx <- cbind2(1, newx)
  #expect_equal(pred_i, newx %*% fit_i$beta[, inds])
});
