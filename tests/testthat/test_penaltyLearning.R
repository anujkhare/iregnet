library(testthat)
library(iregnet)
context("\npenalty learning data set")

data(penalty.learning)

full.fit <- with(penalty.learning, iregnet(
  X.mat, y.mat,
  unreg_sol=FALSE,
  standardize=TRUE,
  maxiter=1e5))

test_that("no lambda=0 when unreg_sol=FALSE", {
  n.lambda.zero <- sum(full.fit$lambda==0)
  expect_equal(n.lambda.zero, 0)
})

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
  lambda=full.fit$lambda,
  unreg_sol=FALSE,
  standardize=TRUE,
  debug=1,
  maxiter=1e5)

test_that("lambda same in fits to full and train data", {
  expect_equal(full.fit$lambda, fit$lambda)
})

test_that("no lambda=0 in fit to train data", {
  n.lambda.zero <- sum(fit$lambda==0)
  expect_equal(n.lambda.zero, 0)
})

if(FALSE){
  conv.txt <- "
iter=1000 lambda=68 beta_19 not converged, abs_change=0.000319 > 0.000100=threshold
iter=1000 lambda=68 beta_22 not converged, abs_change=0.000334 > 0.000100=threshold
iter=1000 lambda=69 beta_19 not converged, abs_change=0.000524 > 0.000100=threshold
iter=1000 lambda=69 beta_22 not converged, abs_change=0.000548 > 0.000100=threshold
iter=1000 lambda=70 beta_19 not converged, abs_change=0.000705 > 0.000100=threshold
iter=1000 lambda=70 beta_22 not converged, abs_change=0.000735 > 0.000100=threshold
iter=1000 lambda=85 beta_9 not converged, abs_change=0.000231 > 0.000100=threshold
iter=1000 lambda=85 beta_22 not converged, abs_change=0.000324 > 0.000100=threshold
iter=1000 lambda=86 beta_3 not converged, abs_change=0.000659 > 0.000100=threshold
iter=1000 lambda=86 beta_8 not converged, abs_change=0.000476 > 0.000100=threshold
iter=1000 lambda=86 beta_9 not converged, abs_change=0.000264 > 0.000100=threshold
iter=1000 lambda=86 beta_11 not converged, abs_change=0.000102 > 0.000100=threshold
iter=1000 lambda=86 beta_13 not converged, abs_change=0.000119 > 0.000100=threshold
iter=1000 lambda=86 beta_18 not converged, abs_change=0.000133 > 0.000100=threshold
iter=1000 lambda=86 beta_22 not converged, abs_change=0.000363 > 0.000100=threshold
iter=1000 lambda=87 beta_3 not converged, abs_change=0.000146 > 0.000100=threshold
iter=1000 lambda=87 beta_8 not converged, abs_change=0.000173 > 0.000100=threshold
iter=1000 lambda=87 beta_9 not converged, abs_change=0.000240 > 0.000100=threshold
iter=1000 lambda=87 beta_22 not converged, abs_change=0.000342 > 0.000100=threshold
iter=1000 lambda=88 beta_0 not converged, abs_change=0.000149 > 0.000100=threshold
iter=1000 lambda=88 beta_1 not converged, abs_change=0.000666 > 0.000100=threshold
iter=1000 lambda=88 beta_3 not converged, abs_change=0.000838 > 0.000100=threshold
iter=1000 lambda=88 beta_8 not converged, abs_change=0.000294 > 0.000100=threshold
iter=1000 lambda=88 beta_9 not converged, abs_change=0.000192 > 0.000100=threshold
iter=1000 lambda=88 beta_11 not converged, abs_change=0.000388 > 0.000100=threshold
iter=1000 lambda=88 beta_13 not converged, abs_change=0.000286 > 0.000100=threshold
iter=1000 lambda=88 beta_14 not converged, abs_change=0.000106 > 0.000100=threshold
iter=1000 lambda=88 beta_18 not converged, abs_change=0.000115 > 0.000100=threshold
iter=1000 lambda=88 beta_22 not converged, abs_change=0.000299 > 0.000100=threshold
iter=1000 lambda=89 beta_0 not converged, abs_change=0.000124 > 0.000100=threshold
iter=1000 lambda=89 beta_1 not converged, abs_change=0.001011 > 0.000100=threshold
iter=1000 lambda=89 beta_3 not converged, abs_change=0.001495 > 0.000100=threshold
iter=1000 lambda=89 beta_5 not converged, abs_change=0.000104 > 0.000100=threshold
iter=1000 lambda=89 beta_8 not converged, abs_change=0.000595 > 0.000100=threshold
iter=1000 lambda=89 beta_9 not converged, abs_change=0.000161 > 0.000100=threshold
iter=1000 lambda=89 beta_11 not converged, abs_change=0.000610 > 0.000100=threshold
iter=1000 lambda=89 beta_13 not converged, abs_change=0.000308 > 0.000100=threshold
iter=1000 lambda=89 beta_14 not converged, abs_change=0.000158 > 0.000100=threshold
iter=1000 lambda=89 beta_15 not converged, abs_change=0.000111 > 0.000100=threshold
iter=1000 lambda=89 beta_18 not converged, abs_change=0.000138 > 0.000100=threshold
iter=1000 lambda=89 beta_20 not converged, abs_change=0.000246 > 0.000100=threshold
iter=1000 lambda=89 beta_22 not converged, abs_change=0.000180 > 0.000100=threshold
iter=1000 lambda=90 beta_0 not converged, abs_change=0.000116 > 0.000100=threshold
iter=1000 lambda=90 beta_1 not converged, abs_change=0.000875 > 0.000100=threshold
iter=1000 lambda=90 beta_3 not converged, abs_change=0.001022 > 0.000100=threshold
iter=1000 lambda=90 beta_4 not converged, abs_change=0.000110 > 0.000100=threshold
iter=1000 lambda=90 beta_8 not converged, abs_change=0.000281 > 0.000100=threshold
iter=1000 lambda=90 beta_9 not converged, abs_change=0.000151 > 0.000100=threshold
iter=1000 lambda=90 beta_11 not converged, abs_change=0.000531 > 0.000100=threshold
iter=1000 lambda=90 beta_13 not converged, abs_change=0.000217 > 0.000100=threshold
iter=1000 lambda=90 beta_14 not converged, abs_change=0.000138 > 0.000100=threshold
iter=1000 lambda=90 beta_15 not converged, abs_change=0.000109 > 0.000100=threshold
iter=1000 lambda=90 beta_18 not converged, abs_change=0.000119 > 0.000100=threshold
iter=1000 lambda=90 beta_20 not converged, abs_change=0.000186 > 0.000100=threshold
iter=1000 lambda=90 beta_22 not converged, abs_change=0.000183 > 0.000100=threshold
iter=1000 lambda=91 beta_1 not converged, abs_change=0.000423 > 0.000100=threshold
iter=1000 lambda=91 beta_3 not converged, abs_change=0.000395 > 0.000100=threshold
iter=1000 lambda=91 beta_4 not converged, abs_change=0.000324 > 0.000100=threshold
iter=1000 lambda=91 beta_9 not converged, abs_change=0.000138 > 0.000100=threshold
iter=1000 lambda=91 beta_11 not converged, abs_change=0.000379 > 0.000100=threshold
iter=1000 lambda=91 beta_14 not converged, abs_change=0.000277 > 0.000100=threshold
iter=1000 lambda=91 beta_20 not converged, abs_change=0.000190 > 0.000100=threshold
iter=1000 lambda=91 beta_22 not converged, abs_change=0.000142 > 0.000100=threshold
iter=1000 lambda=92 beta_3 not converged, abs_change=0.000577 > 0.000100=threshold
iter=1000 lambda=92 beta_4 not converged, abs_change=0.000194 > 0.000100=threshold
iter=1000 lambda=92 beta_8 not converged, abs_change=0.000424 > 0.000100=threshold
iter=1000 lambda=92 beta_9 not converged, abs_change=0.000147 > 0.000100=threshold
iter=1000 lambda=92 beta_13 not converged, abs_change=0.000118 > 0.000100=threshold
iter=1000 lambda=92 beta_14 not converged, abs_change=0.000151 > 0.000100=threshold
iter=1000 lambda=92 beta_18 not converged, abs_change=0.000144 > 0.000100=threshold
iter=1000 lambda=92 beta_20 not converged, abs_change=0.000128 > 0.000100=threshold
iter=1000 lambda=92 beta_22 not converged, abs_change=0.000176 > 0.000100=threshold
iter=1000 lambda=93 beta_1 not converged, abs_change=0.000529 > 0.000100=threshold
iter=1000 lambda=93 beta_3 not converged, abs_change=0.000795 > 0.000100=threshold
iter=1000 lambda=93 beta_4 not converged, abs_change=0.000256 > 0.000100=threshold
iter=1000 lambda=93 beta_8 not converged, abs_change=0.000819 > 0.000100=threshold
iter=1000 lambda=93 beta_9 not converged, abs_change=0.000175 > 0.000100=threshold
iter=1000 lambda=93 beta_11 not converged, abs_change=0.000240 > 0.000100=threshold
iter=1000 lambda=93 beta_13 not converged, abs_change=0.000116 > 0.000100=threshold
iter=1000 lambda=93 beta_18 not converged, abs_change=0.000179 > 0.000100=threshold
iter=1000 lambda=93 beta_22 not converged, abs_change=0.000241 > 0.000100=threshold
iter=1000 lambda=94 beta_1 not converged, abs_change=0.000535 > 0.000100=threshold
iter=1000 lambda=94 beta_3 not converged, abs_change=0.000725 > 0.000100=threshold
iter=1000 lambda=94 beta_4 not converged, abs_change=0.000244 > 0.000100=threshold
iter=1000 lambda=94 beta_8 not converged, abs_change=0.000751 > 0.000100=threshold
iter=1000 lambda=94 beta_9 not converged, abs_change=0.000181 > 0.000100=threshold
iter=1000 lambda=94 beta_11 not converged, abs_change=0.000288 > 0.000100=threshold
iter=1000 lambda=94 beta_13 not converged, abs_change=0.000130 > 0.000100=threshold
iter=1000 lambda=94 beta_18 not converged, abs_change=0.000148 > 0.000100=threshold
iter=1000 lambda=94 beta_22 not converged, abs_change=0.000252 > 0.000100=threshold
iter=1000 lambda=95 beta_1 not converged, abs_change=0.000611 > 0.000100=threshold
iter=1000 lambda=95 beta_3 not converged, abs_change=0.000681 > 0.000100=threshold
iter=1000 lambda=95 beta_4 not converged, abs_change=0.000136 > 0.000100=threshold
iter=1000 lambda=95 beta_8 not converged, abs_change=0.000686 > 0.000100=threshold
iter=1000 lambda=95 beta_9 not converged, abs_change=0.000185 > 0.000100=threshold
iter=1000 lambda=95 beta_11 not converged, abs_change=0.000381 > 0.000100=threshold
iter=1000 lambda=95 beta_13 not converged, abs_change=0.000107 > 0.000100=threshold
iter=1000 lambda=95 beta_18 not converged, abs_change=0.000114 > 0.000100=threshold
iter=1000 lambda=95 beta_22 not converged, abs_change=0.000253 > 0.000100=threshold
iter=1000 lambda=96 beta_1 not converged, abs_change=0.000676 > 0.000100=threshold
iter=1000 lambda=96 beta_3 not converged, abs_change=0.000696 > 0.000100=threshold
iter=1000 lambda=96 beta_4 not converged, abs_change=0.000162 > 0.000100=threshold
iter=1000 lambda=96 beta_8 not converged, abs_change=0.000626 > 0.000100=threshold
iter=1000 lambda=96 beta_9 not converged, abs_change=0.000185 > 0.000100=threshold
iter=1000 lambda=96 beta_11 not converged, abs_change=0.000344 > 0.000100=threshold
iter=1000 lambda=96 beta_13 not converged, abs_change=0.000153 > 0.000100=threshold
iter=1000 lambda=96 beta_18 not converged, abs_change=0.000143 > 0.000100=threshold
iter=1000 lambda=96 beta_22 not converged, abs_change=0.000254 > 0.000100=threshold
iter=1000 lambda=97 beta_0 not converged, abs_change=0.000117 > 0.000100=threshold
iter=1000 lambda=97 beta_1 not converged, abs_change=0.000512 > 0.000100=threshold
iter=1000 lambda=97 beta_3 not converged, abs_change=0.000115 > 0.000100=threshold
iter=1000 lambda=97 beta_4 not converged, abs_change=0.001091 > 0.000100=threshold
iter=1000 lambda=97 beta_5 not converged, abs_change=0.000210 > 0.000100=threshold
iter=1000 lambda=97 beta_8 not converged, abs_change=0.000363 > 0.000100=threshold
iter=1000 lambda=97 beta_9 not converged, abs_change=0.000151 > 0.000100=threshold
iter=1000 lambda=97 beta_10 not converged, abs_change=0.001366 > 0.000100=threshold
iter=1000 lambda=97 beta_11 not converged, abs_change=0.000152 > 0.000100=threshold
iter=1000 lambda=97 beta_14 not converged, abs_change=0.000215 > 0.000100=threshold
iter=1000 lambda=97 beta_15 not converged, abs_change=0.000115 > 0.000100=threshold
iter=1000 lambda=97 beta_20 not converged, abs_change=0.000261 > 0.000100=threshold
iter=1000 lambda=97 beta_22 not converged, abs_change=0.000187 > 0.000100=threshold
iter=1000 lambda=98 beta_3 not converged, abs_change=0.000278 > 0.000100=threshold
iter=1000 lambda=98 beta_4 not converged, abs_change=0.002066 > 0.000100=threshold
iter=1000 lambda=98 beta_5 not converged, abs_change=0.000441 > 0.000100=threshold
iter=1000 lambda=98 beta_8 not converged, abs_change=0.000251 > 0.000100=threshold
iter=1000 lambda=98 beta_10 not converged, abs_change=0.002604 > 0.000100=threshold
iter=1000 lambda=98 beta_12 not converged, abs_change=0.000131 > 0.000100=threshold
iter=1000 lambda=98 beta_13 not converged, abs_change=0.000174 > 0.000100=threshold
iter=1000 lambda=98 beta_14 not converged, abs_change=0.000344 > 0.000100=threshold
iter=1000 lambda=98 beta_15 not converged, abs_change=0.000139 > 0.000100=threshold
iter=1000 lambda=98 beta_20 not converged, abs_change=0.000442 > 0.000100=threshold
iter=1000 lambda=99 beta_0 not converged, abs_change=0.000120 > 0.000100=threshold
iter=1000 lambda=99 beta_1 not converged, abs_change=0.000259 > 0.000100=threshold
iter=1000 lambda=99 beta_3 not converged, abs_change=0.000232 > 0.000100=threshold
iter=1000 lambda=99 beta_4 not converged, abs_change=0.002429 > 0.000100=threshold
iter=1000 lambda=99 beta_5 not converged, abs_change=0.000467 > 0.000100=threshold
iter=1000 lambda=99 beta_8 not converged, abs_change=0.000385 > 0.000100=threshold
iter=1000 lambda=99 beta_10 not converged, abs_change=0.003051 > 0.000100=threshold
iter=1000 lambda=99 beta_13 not converged, abs_change=0.000111 > 0.000100=threshold
iter=1000 lambda=99 beta_14 not converged, abs_change=0.000417 > 0.000100=threshold
iter=1000 lambda=99 beta_15 not converged, abs_change=0.000175 > 0.000100=threshold
iter=1000 lambda=99 beta_20 not converged, abs_change=0.000531 > 0.000100=threshold
iter=1000 lambda=100 beta_0 not converged, abs_change=0.000189 > 0.000100=threshold
iter=1000 lambda=100 beta_1 not converged, abs_change=0.009196 > 0.000100=threshold
iter=1000 lambda=100 beta_2 not converged, abs_change=0.014089 > 0.000100=threshold
iter=1000 lambda=100 beta_3 not converged, abs_change=0.013880 > 0.000100=threshold
iter=1000 lambda=100 beta_4 not converged, abs_change=0.021450 > 0.000100=threshold
iter=1000 lambda=100 beta_5 not converged, abs_change=0.001476 > 0.000100=threshold
iter=1000 lambda=100 beta_6 not converged, abs_change=0.015418 > 0.000100=threshold
iter=1000 lambda=100 beta_7 not converged, abs_change=0.004489 > 0.000100=threshold
iter=1000 lambda=100 beta_8 not converged, abs_change=0.006411 > 0.000100=threshold
iter=1000 lambda=100 beta_9 not converged, abs_change=0.001652 > 0.000100=threshold
iter=1000 lambda=100 beta_10 not converged, abs_change=0.020020 > 0.000100=threshold
iter=1000 lambda=100 beta_11 not converged, abs_change=0.004888 > 0.000100=threshold
iter=1000 lambda=100 beta_12 not converged, abs_change=0.001375 > 0.000100=threshold
iter=1000 lambda=100 beta_13 not converged, abs_change=0.000826 > 0.000100=threshold
iter=1000 lambda=100 beta_14 not converged, abs_change=0.004221 > 0.000100=threshold
iter=1000 lambda=100 beta_15 not converged, abs_change=0.001205 > 0.000100=threshold
iter=1000 lambda=100 beta_16 not converged, abs_change=0.011938 > 0.000100=threshold
iter=1000 lambda=100 beta_17 not converged, abs_change=0.001000 > 0.000100=threshold
iter=1000 lambda=100 beta_18 not converged, abs_change=0.000239 > 0.000100=threshold
iter=1000 lambda=100 beta_19 not converged, abs_change=0.001655 > 0.000100=threshold
iter=1000 lambda=100 beta_20 not converged, abs_change=0.003510 > 0.000100=threshold
iter=1000 lambda=100 beta_21 not converged, abs_change=0.001055 > 0.000100=threshold
iter=1000 lambda=100 beta_22 not converged, abs_change=0.000498 > 0.000100=threshold
"
  library(namedCapture)
  pattern <- paste0(
    "lambda=",
    "(?<lambda>[0-9]+)",
    " ",
    "(?<variable>[^ ]+)",
    " not converged, abs_change=",
    "(?<abs_change>[^ ]+)")
  conv.df <- str_match_all_named(conv.txt, pattern, list(lambda=as.integer, abs_change=as.numeric))[[1]]
  library(Matrix)
  Matrix(with(conv.df, table(variable, lambda)))
}

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
  unreg_sol=FALSE,
  standardize=FALSE,
  debug=1,
  maxiter=1e5)
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

