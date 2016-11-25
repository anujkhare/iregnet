library(testthat)
context("all open intervals\n")
library(iregnet)

set.seed(1)
toy.features <- matrix(rnorm(10), 5, 2)
toy.targets <- rbind(
  c(1, 2),
  c(-Inf, 5),
  c(1, Inf),
  c(-Inf, 3),
  c(-Inf, 4))

test_that("iregnet works with at least one closed interval", {
  fit <- iregnet(toy.features, toy.targets)
  expect_true(inherits(fit, "iregnet"))
})

test_that("iregnet works for all open intervals", {
  fit <- iregnet(toy.features[-1,], toy.targets[-1,])
  expect_true(inherits(fit, "iregnet"))
})

