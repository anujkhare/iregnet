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
    return(list(coef=coef, scale=scale, lambda=fit$lambda, x=x, y=y, type="relax.lasso", family=family))
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
    return(list(coef=ada.beta, scale=fit_lasso$scale, lambda=fit_lasso$lambda, omega=omega, x=x, y=y, type="ada.lasso", family=family))
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

cv.lasso <- function(fit, nfolds=10, foldid=NULL,  rep=1, ...){
    x <- fit$x
    y <- fit$y
    n <- nrow(x)
    num_lambda <- length(fit$lambda)
    if(is.null(foldid)){
        foldid<-sample(rep(1:nfolds, each=ceiling(n/nfolds)), size=n, repl=FALSE) # may happen that not equally distributed
    } else nfolds <- length(unique(foldid))
    loglik <- matrix(0, nrow=nfolds*rep, ncol=num_lambda)
    for(irep in 1:rep){
        for(ifold in 1:nfolds){
            is.train <- foldid!=ifold
            is.valid <- foldid==ifold
            if(fit$type=="lasso") fit.train <- lasso(x=x[is.train, ], y=y[is.train,], ...)
            if(fit$type=="relax.lasso") fit.train <- relax.lasso(x=x[is.train, ], y=y[is.train,], ...)
            if(fit$type=="ada.lasso") fit.train <- ada.lasso(x=x[is.train, ], y=y[is.train,], omega=omega,  ...)
            loglik[ifold,] <- unlist(LAPPLY(1:num_lambda, function(a) 
                survreg(y[is.valid,]~x[is.valid,], dist=fit$family, init=fit.train$coef[,a], 
                        scale=fit.train$scale[a], control=survreg.control(maxiter=0))$loglik[2]))
        }
    }
    cvm <- apply(loglik,2,mean)*nfolds
    cvm.sd <- apply(loglik,2,sd)
    s <- which.max(cvm)
    return(list(lasso.fit = fit, cvm=cvm, cvm.sd=cvm.sd, s=s, lambda.min=fit$lambda[s], type=fit$type, family=fit$family))
}

plot.cv.lasso <- function(fit, adj=5){
    num_lambda <- length(fit$lasso.fit$lambda)
    active <- apply(fit$lasso.fit$coef[-1,], 2, function(X) sum(X!=0))
    xv <- log(fit$lasso.fit$lambda+fit$lasso.fit$lambda[num_lambda-1]*adj)
    ats <- c(xv[1], seq(xv[2], xv[num_lambda], length.out=5)) 
    iats <- c(1, seq(2, num_lambda, length.out=10))
    labs <- round(exp(ats)-fit$lasso.fit$lambda[num_lambda-1]*adj,4)
    ylim <- c(min(fit$cvm - fit$cvm.sd), max(fit$cvm + fit$cvm.sd))
    plot(xv, fit$cvm, type="o", ylim=ylim, axes=F, ylab="CV log Likelihood", xlab="lambda", col="red")
    axis(2)
    for(i in 1:length(fit$lasso.fit$lambda)) lines(c(xv[i], xv[i]), c(fit$cvm[i]-fit$cvm.sd[i], fit$cvm[i]+fit$cvm.sd[i]), col="black") 
    axis(1, at=ats, labels=labs)
    abline(h=fit$cvm[fit$s], lty=3)
    abline(v=xv[fit$s], lty=3)
    points(xv[fit$s], fit$cvm[fit$s], pch=18, col="red")
    axis(3, at=xv[iats], labels=active[iats])
    box()
}

coef.cv.lasso <- function(fit){
    return(fit$lasso.fit$coef[,fit$s])
}

scale.cv.lasso <- function(fit){
    return(fit$lasso.fit$scale[fit$s])
}

predict.cv.lasso <- function(fit, newdata, type="response", ...){
    # uses predict.survreg
    if(missing(newdata)) newdata <- data.frame(fit$lasso.fit$x)
    active <- (coef.cv.lasso(fit)[-1]!=0)
    y <- fit$lasso.fit$y
    x <- data.frame(fit$lasso.fit$x)
    if(sum(active)==0) fit0 <- survreg(fit$lasso.fit$y~1, dist=fit$family, scale=fit$lasso.fit$scale)
    else {
        fit0 <- survreg(as.formula(paste("y~",paste(names(active)[active], collapse="+"))), 
                        data=x,dist=family, init=coef.cv.lasso(fit)[c(TRUE,active)], scale=scale.cv.lasso(fit), control=survreg.control(maxiter=0))
    }
    p<-predict(fit0, newdata, type, ...)
    return(p)
}

surv.cv.lasso <- function(fit, newdata, time, grid=seq(0.001, 0.999, 0.001), se=TRUE){
    quant<-predict.cv.lasso(fit, newdata=newdata, type="quantile", p=grid, se=TRUE)
    return(unlist(sapply(1:nrow(quant$fit), function(X) 1-max(tail(grid[quant$fit[X,]<=time],1),0))))
}

