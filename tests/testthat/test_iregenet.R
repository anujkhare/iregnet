context("iregnet prediction")
library(iregnet)



test_that( "iregnet function converged", {
  load("data/H3K36me3_TDH_immune.Rdata")
  fit <- with(H3K36me3_TDH_immune , iregnet(inputs , outputs))
  expect_equal(fit$error_status, 0)
})


