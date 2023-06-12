lasso <- function(x, y, family,...){
    fit <- iregnet(x=x, y=y, family=family,...)
    return(list(coef=fit$beta, scale=fit$scale, lambda=fit$lambda, x=x, y=y, type="lasso", family=family))
}

relax.lasso <- function(x, y, family, ...){
    require(survival)
    fit <- lasso(x=x, y=y, family=family, ...)
    num_lambda <- length(fit$lambda)
    relaxfit <- LAPPLY(1:num_lambda, function(a) {
        active <- fit$coef[-1,a]!=0
        if(sum(active)==0) ret_fit <- survreg(y~1, dist=family)
        else {
            x_active <- data.frame(x[,active, drop=FALSE])
            ret_fit <- survreg(as.formula(paste("y~",paste(names(active)[active], collapse="+"))), data=x_active,dist=family)
        }
        ret_fit
        })
    coef=matrix(unlist(sapply(1:num_lambda, function(a) {
        beta_tmp<-rep(0, ncol(x)+1)
        names(beta_tmp)<-c("(Intercept)", colnames(x))
        beta_tmp[names(coef(relaxfit[[a]]))]<-coef(relaxfit[[a]])
        beta_tmp
        })), nrow=ncol(x)+1, ncol=num_lambda, byrow=FALSE)
    rownames(coef) <- c("(Intercept)", colnames(x))
    scale=sapply(1:num_lambda, function(a) {
            relaxfit[[a]]$scale
        })
    return(list(coef=coef, scale=scale, lambda=fit$lambda, x=x, y=y, type="relax.lasso"))
    }


ada.lasso <- function(x, y, family, omega=1, ...){
    require(survival)
    fit_init <- survreg(y~x, dist=family)
    weights <- abs(coef(fit_init)[-1])**omega
    xw <- x %*% diag(weights)
    fit_lasso <- iregnet(y=y, x=xw, standardize=FALSE, family=family, ...)
    ada.beta<-fit_lasso$beta
    ada.beta <- rbind(ada.beta[1,], diag(weights) %*% fit_lasso$beta[-1,])
    rownames(ada.beta) <- c("(Intercept)",colnames(x))
    return(list(coef=ada.beta, scale=fit_lasso$scale, lambda=fit_lasso$lambda, omega=omega, x=x, y=y, type="ada.lasso"))
}

plotlasso <- function(fit, xvar="log"){
    if(xvar=="log") {
        xv <- log(fit$lambda+1)
        xlabel <- "log(lambda+1)"
    } else if(xvar=="L1norm") {
        xv <- apply(fit$coef[-1,],2,function(X) sum(abs(X)))
        xlabel <- "L1 norm"
    }
    ylim <- range(fit$coef)
    xlim <- range(xv)
    plot(x=xv, y=fit$coef[1,], type="l", lty=3, xlim=xlim, ylim=ylim, ylab="beta", xlab=xlabel)
    for(j in 1:nrow(fit$coef[-1,])) lines(x=xv, y=fit$coef[j+1,], type="l", lty=1, col=j)
}

cv.lasso <- function(fit, nfolds=10, foldid=NULL,  ...){
    x <- fit$x
    y <- fit$y
    n <- nrow(x)
    num_lambda <- length(fit$lambda)
    if(is.null(foldid)){
        foldid<-sample(rep(1:nfolds, each=ceiling(n/nfolds)), size=n, repl=FALSE) # may happen that not equally distributed
    } else nfolds <- length(unique(foldid))
    loglik <- matrix(0, nrow=nfolds, ncol=num_lambda)
    for(ifold in 1:nfolds){
        is.train <- foldid!=ifold
        is.valid <- foldid==ifold
        if(fit$type=="lasso") fit.train <- lasso(x=x[is.train, ], y=y[is.train,], family=fit$family, ...)
        if(fit$type=="relax.lasso") fit.train <- relax.lasso(x=x[is.train, ], y=y[is.train,], family=fit$family, ...)
        if(fit$type=="ada.lasso") fit.train <- ada.lasso(x=x[is.train, ], y=y[is.train,], omega=omega, family=fit$family, ...)
        loglik[ifold,] <- LAPPLY(1:num_lambda, function(a) survreg(y[is.valid,]~x[is.valid], dist=family, init=fit.train$coef[,a], scale=fit.train$scale[a], estimate.scale=FALSE, control=survreg.control(maxiter=0)))
    }
    cvm <- apply(loglik,2,sum)
    cvm.sd <- apply(loglik,2,sd)
    s <- fit$lambda[which.max(cvm)]
    return(list(lasso.fit = fit, cvm=cvm, cvm.sd=svm.sd, s=s, lambda.min=fit$lambda[s], type=fit$type, family=fit$family))
}
