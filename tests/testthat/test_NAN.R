library(testthat)
context("NAN")
library(iregnet)
data(realNAN)
set.seed(1)
fit <- cv.iregnet(
  realNAN$feature.mat, realNAN$target.mat,
  family="gaussian")
