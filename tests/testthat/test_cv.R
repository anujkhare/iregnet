library(testthat)
context("cv\n")
library(iregnet)
library(data.table)
data(penalty.learning)
fit <- with(penalty.learning, cv.iregnet(X.mat, y.mat, family="gaussian"))

test_that("plot(cv.iregnet result) yields ggplot", {
  result <- plot(fit)
  expect_true(is.ggplot(result))
})

