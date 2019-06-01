library(testthat)
context("\nElemStatsLearn dataset")
library(iregnet)
library(survival)
library(glmnet)

test_that("ElemStatsLearn data - coefficients are calculated correctly wrt survival and glmnet:", {
  alpha <- 0.6

  data(prostate,package="ElemStatLearn")
  pros <- subset(prostate,select=-train,train==TRUE)
  ycol <- which(names(pros)=="lpsa")
  X.unscaled <- as.matrix(pros[-ycol])
  y.unscaled <- pros[[ycol]]
  M <- matrix(
    colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
  X.centered <- X.unscaled - M
  sd.vec <- apply(X.unscaled, 2, sd)
  S <- diag(1/sd.vec)
  X.scaled <- X.centered %*% S
  dimnames(X.scaled) <- dimnames(X.unscaled)
  m <- mean(y.unscaled)
  sigma <- sd(y.unscaled)
  y.scaled <- (y.unscaled - m)/sigma

  X <- X.scaled
  y <- y.scaled

  fit_s <- survreg(Surv(y, rep(1, length(y))) ~ X, dist = "gaussian")
  fit_i <- iregnet(X, cbind(y, y), "gaussian", maxiter=1e5, thresh=1e-7, standardize=F, alpha=alpha, scale=1, estimate_scale=F)

  lambda_path <- fit_i$lambda * (fit_i$scale ** 2)

  fit_g <- glmnet(X, y, "gaussian", lambda = lambda_path, standardize=F, maxit=1e5, thresh=1e-7, alpha=alpha)

  expect_equal(as.double(fit_s$coefficients), as.double(fit_i$beta[, fit_i$num_lambda]), tolerance = 1e-3)
  expect_equal(as.double(fit_i$beta), as.double(coef(fit_g)), tolerance=1e-3)
})
