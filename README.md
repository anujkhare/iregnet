[![Travis-CI Build Status](https://travis-ci.org/anujkhare/iregnet.svg?branch=master)](https://travis-ci.org/anujkhare/iregnet)

This is an R package, in development, for regularized interval regression.

## Installation
To install this version of the package, use the devtools package:
```R
devtools::install_github("georgheinze/iregnet")
```

## Example usage:
```R
# y can be a 2 column matrix.
X <- matrix(rnorm(50), 10, 5)
y <- matrix(rnorm(20), 10, 2)
y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
fit1 <- iregnet(X, y)

# Surv objects from survival are also supported.
data("ovarian")
X <- cbind(ovarian$ecog.ps, ovarian$rx)
y <- Surv(ovarian$futime, ovarian$fustat)
fit2 <- iregnet(X, y)

# Log-Gaussian is same as Gaussian with log-transformed data
X <- matrix(rnorm(50), 10, 5)
y <- matrix(abs(rnorm(20)), 10, 2)
y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
fit3 <- iregnet(X, log(y), "gaussian")
fit4 <- iregnet(X, y, "loggaussian")

# Scale parameter can be fixed by setting the estimate_scale flag.
X <- matrix(rnorm(50), 10, 5)
y <- matrix(rnorm(20), 10, 2)
y <- t(apply(y, 1, sort)) # intervals must be non-decreasing
fit5 <- iregnet(X, y, scale_init=1, estimate_scale=F)
```

Detailed documentation of each parameter is provided in R help for iregnet.

For more details about the package, visit the following link:
https://github.com/rstats-gsoc/gsoc2016/wiki/Regularized-interval-regression
