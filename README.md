[![Travis-CI Build Status](https://travis-ci.org/anujkhare/iregnet.svg?branch=master)](https://travis-ci.org/anujkhare/iregnet)

This is a R package, in development, for regularized interval regression.

To install the package, use the devtools package:
```R
devtools::install_github("anujkhare/iregnet")
```

Example usage:
```R
library(iregnet)
data(ovarian)

# y should be two-column, with NA denoting censoring
y_l <- y_r <- ovarian$futime
y_r[ovarian$fustat == 0] <- NA
y <- cbind(y_l, y_r)
X <- cbind(ovarian$ecog.ps, ovarian$rx)

fit <- iregnet(X, y, family="gaussian", alpha=1, intercept = T)
print(fit)
```

For more details, visit the following link:
https://github.com/rstats-gsoc/gsoc2016/wiki/Regularized-interval-regression
sdfh
