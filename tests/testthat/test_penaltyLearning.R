library(testthat)
library(iregnet)
context("\npenalty learning data set")

data(penalty.learning)

chrom.vec <- sub(":.*", "", rownames(penalty.learning$X.mat))
table(chrom.vec)
train.chroms <- c("chr1", "chr9")
sets <-
  list(train=chrom.vec %in% train.chroms,
       validation=! chrom.vec %in% train.chroms)
X.train <- penalty.learning$X.mat[sets$train,]
y.train <- penalty.learning$y.mat[sets$train,]
fit <- iregnet(
  X.train, y.train,
  ##threshold=1e-1,
  standardize=TRUE,
  debug=1)

test_that("predict function same as matrix multiplication when standardize=TRUE", {
  expect_equal(cbind(1, X.train) %*% fit$beta, predict(fit, X.train))
})

l1.norm.vec <- colSums(abs(fit$beta[-1,]))
test_that("only one model should have L1 arclength 0", {
  n.zero <- sum(l1.norm.vec == 0)
  expect_equal(n.zero, 1)
})

pred.loglik <- function(pred.mean, pred.scale, target.mat){
  stopifnot(identical(ncol(target.mat), 2L))
  n.obs <- nrow(target.mat)
  stopifnot(identical(nrow(pred.mean), n.obs))
  stopifnot(identical(nrow(pred.scale), n.obs))
  prob.below.left <- pnorm(target.mat[,1], pred.mean, pred.scale)
  prob.below.right <- pnorm(target.mat[,2], pred.mean, pred.scale)
  log(prob.below.right-prob.below.left)
}

set.error.list <- list()
for(set.name in names(sets)){
  is.set <- sets[[set.name]]
  pred.mat <- predict(fit, penalty.learning$X.mat[is.set,])
  y.set <- penalty.learning$y.mat[is.set,]
  too.lo <- colSums(pred.mat < y.set[,1])
  too.hi <- colSums(y.set[,2] < pred.mat)
  errors <- too.lo + too.hi
  scale.mat <- matrix(fit$scale, nrow(y.set), length(fit$scale), byrow=TRUE)
  log.prob <- pred.loglik(pred.mat, scale.mat, y.set)
  set.error.list[[set.name]] <- data.frame(
    set.name,
    lambda=fit$lambda,
    too.lo,
    too.hi,
    errors,
    iterations=fit$n_iters,
    surrogate.loss=-colMeans(log.prob),
    percent.error=100*errors/nrow(y.set))
}
set.error <- do.call(rbind, set.error.list)
train.error <- subset(set.error, set.name=="train")

if(interactive()){
  library(ggplot2)

  tidy <- tidydf(fit)
  some <- subset(tidy, variable != "(Intercept)" & arclength < 0.25)

  ggplot()+
    geom_path(aes(-log(lambda), weight, color=variable, group=variable),
              data=some)

  ggplot()+
    geom_path(aes(arclength, weight, color=variable, group=variable),
              data=some)+
    geom_point(aes(arclength, weight, color=variable, group=variable),
               shape=1,
               data=subset(some, weight != 0))

  ggplot()+
    geom_line(aes(arclength, weight, color=variable, group=variable),
              data=some)+
    geom_point(aes(arclength, weight, color=variable, group=variable),
               shape=1,
               data=subset(some, weight != 0))

  ggplot()+
    ggtitle("iregnet on penalty.learning data set")+
    ylab("")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(metric ~ ., scales="free")+
    geom_path(aes(-log(lambda), surrogate.loss, color=set.name, group=set.name),
              data=data.frame(set.error, metric="surrogate loss"))+
    geom_path(aes(-log(lambda), iterations),
              data=data.frame(train.error, metric="iterations"))+
    geom_path(aes(-log(lambda), percent.error, color=set.name, group=set.name),
              data=data.frame(set.error, metric="percent incorrect intervals"))

  plot(log(fit$lambda))

}

test_that("train set surrogate loss always decreases", {
  diff.vec <- diff(train.error$surrogate.loss)
  expect_true(all(diff.vec < 0))
})

test_that("error_status is 0", {
  expect_equal(fit$error_status, 0)
})

## Another fit with standardize=FALSE.
X.unscaled <- X.train
mean.vec <- colMeans(X.unscaled)
M <- matrix(
  mean.vec, nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
X.centered <- X.unscaled - M
sd.vec <- sqrt(colSums(X.centered * X.centered)/nrow(X.centered))
S <- diag(1/sd.vec)
X.scaled <- X.centered %*% S
dimnames(X.scaled) <- dimnames(X.unscaled)
ufit <- iregnet(
  X.scaled, y.train,
  ## scale_init=1, estimate_scale=FALSE,
  standardize=FALSE)
pred.mat <- predict(ufit, X.scaled)
test_that("predict function same as matrix multiplication when standardize=FALSE", {
  expect_equal(cbind(1, X.scaled) %*% ufit$beta, pred.mat)
})
## If x is an unscaled p-vector then z = S'(x-m) is a scaled p-vector,
## where m is a p-vector of means and S is a (p x p) diagonal matrix
## of 1/sd values. The iregnet(standardize=FALSE) function gives us a
## p-vector of weights w and an intercept b, so the final prediction
## function is w'z + b = w'S(x-m) + b = w'Sz + b - w'Sm. Below we
## compute the scaled.beta.mat which contains the w'S vectors, and the
## pred.weight.mat which can be used for prediction with unscaled
## data.
scaled.beta.mat <- S %*% ufit$beta[-1,]
intercept.vec <- ufit$beta[1,] - mean.vec %*% scaled.beta.mat
pred.weight.mat <- rbind(intercept.vec, scaled.beta.mat)
my.pred.mat <- cbind(1, X.train) %*% pred.weight.mat
test_that("lambda seq is the same", {
  expect_equal(fit$lambda, ufit$lambda)
})
lambda.diff <- diff(fit$lambda)
test_that("lambda seq decreases", {
  expect_true(all(lambda.diff < 0))
})
test_that("learned weights are the same", {
  expect_equal(unname(pred.weight.mat), unname(fit$beta))
})

