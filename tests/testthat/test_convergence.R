context("iregnet prediction")
library(iregnet)


test_that("iregnet function converged", {
  
  data( H3K36me3_TDH_immune , package = "iregnet")
  fit <- with(H3K36me3_TDH_immune , iregnet::iregnet(inputs , outputs , debug = 1))
  expect_equal(fit$error_status, 0)
})


