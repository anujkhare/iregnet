context("iregnet prediction")
library(iregnet)


## estimated scale - does not converg
test_that("iregnet function converged", {

  data( H3K36me3_TDH_immune , package= "iregnet")
  fit_est <- with(H3K36me3_TDH_immune , iregnet::iregnet(inputs , outputs , debug = 1))
  expect_equal(fit$error_status, 0)
})


## Fixed scale
test_that("iregnet function converged", {
  
  data( H3K36me3_TDH_immune , package= "iregnet")
  fit_fixed <- with(H3K36me3_TDH_immune , iregnet::iregnet(inputs , outputs ,scale_init= 1 ,
                                                     estimate_scale=FALSE ,debug = 1))
  expect_equal(fit$error_status, 0)
})
