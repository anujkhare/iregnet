library(testthat)
context("NAN")
library(iregnet)
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
