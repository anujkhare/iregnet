library(testthat)
context("NAN")
library(iregnet)
require(data.table)
data("neuroblastomaProcessed")
data(realNAN)
set.seed(1)

## No problem here:
fit <- cv.iregnet(
  realNAN$feature.mat[, -27], realNAN$target.mat,
  family="gaussian")

## Error: NANs produced
fit <- cv.iregnet(
  realNAN$feature.mat, realNAN$target.mat,
  family="gaussian")


#Custom data which works fine
test_that("Produces failed to converge warning",{
  set.seed(125)
  
  x <- matrix(rnorm(25), 5, 5)
  y <- rbind(
    c(5, 10),
    c(6, Inf),
    c(-Inf, 2),
    c(1, 2),
    c(-Inf, 3))
  expect_warning(iregnet(x, y), "Failed to converge. Try again after adding more data.")
  
  y <- rbind(
    c(5, 10),
    c(6, Inf),
    c(9, Inf),
    c(3, 4),
    c(10, Inf))
  expect_warning(iregnet(x, y), "Failed to converge. Try again after adding more data.")

  y <- rbind(
    c(1, 2),
    c(-Inf, 5),
    c(3, Inf),
    c(-Inf, 3),
    c(-Inf, 4))
  expect_warning(iregnet(x, y), "Failed to converge. Try again after adding more data.")
});


test_that("Check left and right censorship", {
  x <- matrix(rnorm(25), 5, 5)
  y <- rbind(  
    c(-Inf, 2),
    c(-Inf, 3),
    c(-Inf, 4),
    c(-Inf, 5),
    c(-Inf, 8))
  expect_error(iregnet(x, y), "Target matrix completely left censored. Try adding more data")
  expect_error(iregnet(x, Surv(y[,1], y[,2], type = "interval2")), "Target Surv object completely left censored. Try adding more data")
  
  y <- rbind(
    c(3, Inf),
    c(4, Inf),
    c(6, Inf),
    c(1, Inf),
    c(2, Inf))
  expect_error(iregnet(x, y), "Target matrix completely right censored. Try adding more data")
  expect_error(iregnet(x, Surv(y[,1], y[,2], type = "interval2")), "Target Surv object completely right censored. Try adding more data")
});

test_that("NA's and Inf, all the same", {
  set.seed(125)
  x <- matrix(rnorm(25), 5, 5)
  y <- rbind(
    c(5, 10),
    c(6, NA),
    c(-Inf, 2),
    c(1, NA),
    c(-Inf, 3))
  expect_warning(iregnet(x, y), "Failed to converge. Try again after adding more data.")
})
