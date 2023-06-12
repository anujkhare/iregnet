distribution.fun.suffixes <- c(
  gaussian="norm",
  logistic="logis",
  exponential="exponential",
  weibull="weibull")
dexponential <- function(x, m, s, log){
  dexp(x, m, log)
}
pexponential <- function(x, m, s){
  pexp(x, m)
}


paste0get <- function(suffix, prefix){
  get(paste0(prefix, suffix))
}
pfun.list <- lapply(distribution.fun.suffixes, paste0get, "p")
dfun.list <- lapply(distribution.fun.suffixes, paste0get, "d")

##' Compute log-likelihood of iregnet predictions, used to compute the
##' surrogate loss on the validation set in cross-validation
##' (cv.iregnet).
##' Note that for family="weibull" it will compute the likelihood only as sum, not as components
##' @title Compute log likelihood
##' @param y.mat numeric matrix of output variables (n x 2)
##' @param pred.mean numeric matrix of mean parameters (predicted
##'   values, n x p)
##' @param pred.scale numeric matrix of scale parameters (n x p)
##' @param family gaussian, logistic, weibull
##' @export
##' @return numeric matrix of log-likelihood values (n x p)
##' @author Georg Heinze
compute.loglik.weibull <- function(y.mat, pred.mean, pred.scale, family="weibull"){
    compute.loglik(y.mat=y.mat, pred.mean=pred.mean, pred.scale=pred.scale, family="weibull")
}
    

##' Compute log-likelihood of iregnet predictions, used to compute the
##' surrogate loss on the validation set in cross-validation
##' (cv.iregnet).
##' Note that for family="weibull" it will compute the likelihood only as sum, not as components
##' @title Compute log likelihood
##' @param y.mat numeric matrix of output variables (n x 2)
##' @param pred.mean numeric matrix of mean parameters (predicted
##'   values, n x p)
##' @param pred.scale numeric matrix of scale parameters (n x p)
##' @param family gaussian, logistic, weibull
##' @export
##' @return numeric matrix of log-likelihood values (n x p)
##' @author Toby Dylan Hocking
compute.loglik <- function(y.mat, pred.mean, pred.scale, family){
    require(survival)
  stopifnot_error(
    "y.mat must be a 2-column numeric matrix",
    is.matrix(y.mat),
    is.numeric(y.mat),
    ncol(y.mat)==2L)
  n.obs <- nrow(y.mat)
  stopifnot_error(
    "pred.mean must be a numeric matrix with same number of rows as y.mat",
    is.matrix(pred.mean),
    is.numeric(pred.mean),
    nrow(pred.mean)==n.obs)
  stopifnot_error(
    "pred.scale must be a numeric matrix with same number of rows as y.mat",
    is.matrix(pred.scale),
    is.numeric(pred.scale),
    nrow(pred.scale)==n.obs)
  stopifnot_error(
    paste("family must be one of", paste(names(pfun.list), collapse=", ")),
    family %in% names(pfun.list))
  if(family=="weibull"){
      status <- (y.mat[,1]==y.mat[,2])
      loglik.mat <- sapply(1:ncol(pred.mean), function(ilambda) 
          survreg(Surv(y.mat[,1],status)~pred.mean[,ilambda]-1, scale=pred.scale[ilambda], init=1, control=survreg.control(maxiter=0))$loglik[2])
  } else {
      pfun <- pfun.list[[family]]
      dfun <- dfun.list[[family]]
      use.density <- y.mat[,1]==y.mat[,2]
      loglik.mat <- matrix(NA, n.obs, ncol(pred.scale))
      loglik.mat[use.density,] <- dfun(
        y.mat[use.density, 1],
        pred.mean[use.density, ],
        pred.scale[use.density, ],
        log=TRUE)
      loglik.mat[!use.density,] <- log(
        pfun(y.mat[!use.density, 2],
             pred.mean[!use.density, ],
             pred.scale[!use.density, ])
        -pfun(y.mat[!use.density, 1],
              pred.mean[!use.density, ],
              pred.scale[!use.density, ]))
  }
  loglik.mat
}

##' First fit iregnet to the entire data set, then use K-fold
##' cross-validation to estimate the validation error of each lambda
##' parameter.
##' @title Cross-validation for iregnet
##' @param x numeric matrix of input features (n x p)
##' @param y numeric matrix of output variables (n x 2)
##' @param family gaussian, logistic, weibull
##' @param nfolds positive integer between 2 and n, by default 10.
##' @param foldid integer vector of length n (fold for each observation), by default we use nfolds.
##' @param ... passed to iregnet for both the full and cv fits.
##' @export
##' @import future.apply
##' @return model fit list of class "cv.iregnet"
##' @author Toby Dylan Hocking, Georg Heinze
cv.iregnet <- function(x, y, family, nfolds, foldid, relax=FALSE, ...){
  if(missing(foldid)){
    if(missing(nfolds)){
      nfolds <- 10L
    }else{
      stopifnot_error(
        "nfolds should be an integer between 2 and nrow(x)",
        is.integer(nfolds),
        2 <= nfolds, nfolds <= nrow(x))
    }
    foldid <- sample(rep(1:nfolds, l=nrow(x)))
  }else{
    stopifnot_error(
      "foldid should be an integer vector of length nrow(x)",
      is.integer(foldid), length(foldid) == nrow(x))
  }

  # If no column names are given
  if(!length(colnames(x))){
    colnames(x) <- paste('x', 1: ncol(x), sep='')
  }

  big.fit <- iregnet(x, y, family=family, unreg_sol=FALSE, ...)
  validation.fold.vec <- unique(foldid)
  
  LAPPLY <- if(requireNamespace("future.apply")){
    future.apply::future_lapply
  }else{
    lapply
  }  

  oldw <- getOption("warn") # Supress warnings for now
  options(warn = -1)
  error <- c() # Vector to record the fold errors

  fold.mat.list <- LAPPLY(validation.fold.vec, function(validation.fold){
    is.validation <- foldid == validation.fold
    is.train <- !is.validation
    x.train <- x[is.train,]
    y.train <- y[is.train,]    
    fit <- iregnet(
      x.train, y.train,
      family=family,
      lambda=big.fit$lambda,
      unreg_sol=FALSE, ...)
    if(relax){
        status<-(y.train[,1]==y.train[,2])*1 ## only right-censoring with relax!
        time <- y.train[,1]
        fitrelax <- LAPPLY(1:length(fit$lambda), function(a) {
                activeset <- fit$beta[-1,a]!=0
#                if(sum(activeset)==0)  formula <- Surv(time,status)~1                 else formula <- as.formula(paste("Surv(time,status)~", paste(colnames(x.train[,activeset]),collapse="+")))
                dataset <- data.frame(time=y.train[,1], status=status, x.train[,activeset,drop=FALSE])
                fit_aft <- survreg(Surv(time,status)~., dataset, dist="weibull")
                fit_aft
            }
            )
        fit$beta <- matrix(unlist(sapply(1:length(fitrelax), function(a) {
            beta_tmp<-rep(0, ncol(x.train)+1)
            names(beta_tmp)<-c("(Intercept)", colnames(x.train))
            beta_tmp[names(coef(fitrelax[[a]]))]<-coef(fitrelax[[a]])
            beta_tmp
        })), nrow=nrow(fit$beta), ncol=length(fitrelax), byrow=FALSE)
        fit$scale <- unlist(sapply(1:length(fitrelax), function(a) fitrelax[[a]]$scale))
        fit$loglik <- unlist(sapply(1:length(fitrelax), function(a) fitrelax[[a]]$loglik[2]))
        rownames(fit$beta) <- names(beta)
    }
    error <<- append(error, fit$error_status)
    pred.center.mat <- predict(fit, x)
    pred.scale.mat <- matrix(
      fit$scale,
      nrow(pred.center.mat),
      ncol(pred.center.mat),
      byrow=TRUE)
    if(family != "weibull") loglik <- compute.loglik(y, pred.center.mat, pred.scale.mat, family=family)
    else {
        loglik.weibull <- cbind(compute.loglik.weibull(y[is.train,], pred.center.mat[is.train,], pred.scale.mat[is.train,], family=family),
                            compute.loglik.weibull(y[is.validation,], pred.center.mat[is.validation,], pred.scale.mat[is.validation,], family=family))
        colnames(loglik.weibull)<-c("train", "validation")
    }
    sets <- list(validation=is.validation, train=is.train)
    one.fold.mat <- matrix(NA, length(fit$lambda), 2, dimnames=list(
      lambda=fit$lambda, set=c("train", "validation")))
    for(set.name in names(sets)){
      is.set <- sets[[set.name]]
      if(family != "weibull") one.fold.mat[, set.name] <- -colMeans(loglik[is.set, ])
      else one.fold.mat[, set.name] <- -loglik.weibull[,set.name]/sum(is.set) # mean loglik instead of sum of components
    }
    one.fold.mat
  })

  options(warn = oldw) # Turn on the warnings
  if(is.element(-1, error))
    warning("Ran out of iterations and failed to converge on some folds. Check the error status in fit.")
  if(is.element(-3, error))
    warning("Failed to converge. Try again after adding more data. Check the error status in fit.")

  fold.array <- array(
    unlist(fold.mat.list),
    c(length(big.fit$lambda), 2, nfolds),
    list(
      lambda=big.fit$lambda,
      set=c("train", "validation"),
      validation.fold=validation.fold.vec))
  mean.mat <- apply(fold.array, c(1,2), mean)
  sd.mat <- apply(fold.array, c(1,2), sd)
  stats.df.list <- list()
  for(set.name in c("train", "validation")){
    stats.df.list[[set.name]] <- data.table::data.table(
      set.name,
      lambda=big.fit$lambda,
      mean=mean.mat[, set.name],
      sd=sd.mat[, set.name])
  }
  stats.df <- do.call(rbind, stats.df.list)
  lik.df.list <- list()
  for(fold.i in seq_along(validation.fold.vec)){
    validation.fold <- validation.fold.vec[[fold.i]]
    one.fold.mat <- fold.mat.list[[fold.i]]
    for(set.name in c("train", "validation")){
      lik.df.list[[paste(fold.i, set.name)]] <- data.table::data.table(
        validation.fold,
        set.name,
        lambda=big.fit$lambda,
        neg.loglik=one.fold.mat[, set.name])
    }
  }
  lik.df <- do.call(rbind, lik.df.list)
  i.min <- which.min(mean.mat[, "validation"])
  min.upper <- mean.mat[i.min, "validation"]+sd.mat[i.min, "validation"]
  is.within <- mean.mat[, "validation"] < min.upper
  i.1sd <- which(is.within)[1]
  big.fit$selected <- c(min=as.integer(i.min), "1sd"=as.integer(i.1sd))
  selected.df <- data.table::data.table(
    set.name="validation",
    type=names(big.fit$selected),
    lambda=big.fit$lambda[big.fit$selected],
    neg.loglik=mean.mat[big.fit$selected, "validation"])
  not.intercept <- big.fit$beta[-1,]
  nonzero.df <- data.table::data.table(
    arclength=colSums(abs(not.intercept)),
    lambda=big.fit$lambda,
    nonzero=colSums(not.intercept != 0))
  weight.and.intercept <- tidydf(big.fit)
  is.intercept <- weight.and.intercept$variable == "(Intercept)"
  weight.df <- weight.and.intercept[!is.intercept, ]
  big.fit$plot.data <- list(
    weights=weight.df,
    stats=stats.df,
    likelihood=lik.df,
    selected=selected.df)
  big.fit$error_status <- error
  class(big.fit) <- c("cv.iregnet", "iregnet", "list")
  big.fit
}