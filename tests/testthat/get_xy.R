get_xy <- function(n_obs, n_vars, type = c('right', 'left', 'none', 'interval'), standardize=std, positive=F) {
  type <- match.arg(type)

  x <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars)
  y <- rnorm(n_obs)
  y_u <- rnorm(n_obs)

  # standardize x and y
  if (standardize == T) {
    for (i in 1:ncol(x)) {
      x[, i] <- (x[, i] - mean(x[, i])) / sd(x[, i]);
    }

    y <- (y - mean(y))
    y <- y / sd(y)

    if (type == 'interval') {
      y_u <- (y_u - mean(y_u))
      y_u <-  y_u / sd(y_u)
    }
  }

  if (positive) {
    y <- abs(y)
    y_u <- abs(y_u)
  }

  if (type == "none") {
    status = rep(1, length(y))
    y_surv <- Surv(time = y, event = status, type = "right")

  } else if (type == 'interval') {
    status <- sample(c(0, 1, 2, 3), size=n_obs, replace=T)
    y <- cbind(y, y_u)
    y <- t(apply(y, 1, sort))    # make sure intervals are increasing in time
    y_surv <- Surv(time=y[, 1], time2=y[, 2], event = status, type = 'interval')

  } else {    # left or right censored
    status <- sample(c(0, 1), size=n_obs, replace=T)
    y_surv <- Surv(time = y, event = status, type = type)
  }

  # get the y matrix
  y <- cbind(y, y)
  if (type=="right") {
    y[status == 0, 2] = NA
  } else if (type == 'interval') {
    y <- NA  # NOTE: Not implemented, not needed, use the Surv object!
  } else {
    y[status == 0, 1] = NA
  }

  return (list("x" = x, "y" = y, "surv" = y_surv))
}
