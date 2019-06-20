library(testthat)
context("all open intervals\n")
library(iregnet)
library(penaltyLearning)

set.seed(1)
toy.features <- matrix(rnorm(10, 1, 0.1), 5, 2)
toy.targets <- rbind(
  c(1, 2),
  c(-Inf, 5),
  c(3, Inf),
  c(-Inf, 3),
  c(-Inf, 4))
if(interactive()){
  colnames(toy.targets) <- c("min.L", "max.L")
  toy.df <- data.frame(toy.targets, toy.features)
  library(survival)
  fit <- survreg(
    Surv(min.L, max.L, type="interval2") ~ X1 + X2,
    toy.df, dist="gaussian")
### This gives a warning but still returns a valid model.
### In survreg.fit(X, Y, weights, offset, init = init, controlvals = control,  :
### Ran out of iterations and did not converge
  pred.vec <- predict(fit)
  is.lo <- pred.vec < toy.targets[,1]
  is.hi <- toy.targets[,2] < pred.vec
  is.error <- is.lo | is.hi
  n.errors <- sum(is.error)
  test_that("survival-2.41.0 predicts correctly for toy data", {
    expect_equal(n.errors, 0)
  })

  test_that("iregnet works with at least one closed interval", {
    fit <- iregnet(toy.features, toy.targets)
    expect_true(inherits(fit, "iregnet"))
  })

  test_that("iregnet works for all open toy intervals", {
    fit <- iregnet(toy.features[-1,], toy.targets[-1,])
    expect_true(inherits(fit, "iregnet"))
  })

  print(.libPaths())

}

if(require(penaltyLearning)){

  data(neuroblastomaProcessed, package="penaltyLearning")

  if(interactive()){
    train.dt <- with(neuroblastomaProcessed, data.table(
      log2.n=feature.mat[, "log2.n"],
      log.hall=feature.mat[, "log.hall"],
      target.mat))
    fit <- survreg(
      Surv(min.L, max.L, type="interval2") ~ log2.n + log.hall,
      train.dt, dist="gaussian")
    ## survival returns a valid model for these data.
  }

  test_that("valid model for two features in neuroblastoma data", {
    fit <- with(neuroblastomaProcessed, iregnet(
      feature.mat[, c("log2.n", "log.hall")],
      target.mat, family="gaussian"))
    expect_is(fit, "iregnet")
    df <- tidydf(fit)
    #ggplot()+
     # geom_point(aes(-log(lambda), weight, color=variable), data=df)
  })

  test_that("valid model for all features in neuroblastoma data", {
    fit <- with(neuroblastomaProcessed, {
      iregnet(feature.mat, target.mat, family="gaussian")
    })
    expect_is(fit, "iregnet")
  })

  sd.vec <- apply(neuroblastomaProcessed$feature.mat, 2, sd)
  is.constant <- sd.vec == 0

  test_that("valid model for non-constant features in neuroblastoma data", {
    fit <- with(neuroblastomaProcessed, {
      iregnet(feature.mat[, !is.constant], target.mat, family="gaussian")
    })
    expect_is(fit, "iregnet")
  })

}
