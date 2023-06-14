#' @title Fit lasso for right-censored AFT model
#' @export
#' @description
#' This is a simple wrapper function to fit accelerated failure time models 
#' using elastic net penalized maximum likeihood, calling \code{iregnet}.
#' Supports gaussian, logistic, Weibull and extreme value distributions.
#'
#' @param x Input matrix of covariates with dimension n_obs * n_vars, with
#' \eqn{nvars \ge 2}. Sparse matrices are not supported.
#' 
#' @param y Response variable. Here, only \code{\link{Surv}} object is supported.
#'
#' @param family The distribution to fit. It can be one of "gaussian", "logistic",
#' "loggaussian", "loglogistic", "extreme_value", "exponential", "weibull". Partial matching
#' is allowed.
#' \cr \emph{Default: "gaussian"}
#'
#' @param ... Further arguments passed to \code{\link{iregnet}} 
#' 
#' 
#' @details This is a wrapper function for \code{\link{iregnet}} with simpler
#' output, facilitating simple cross-validation.
#' @return Returns an object  with the following elements:\cr
#' \tabular{ll}{
#'  \code{coef} \tab Matrix of size \code{(n_vars+1) * num_lambda} containing
#'  intercept, coefficients of \code{X} for each \code{lambda} in the fit model.
#'    \cr
#'  \code{scale} \tab Vector of size \code{num_lambda} of estimated
#'    scale at each \code{lambda} value. \cr
#'  \code{lambda} \tab Vector of size \code{num_lambda} of (calculated or
#'   supplied) regularization parameter \code{lambda} values. \cr
#'  \code{x} \tab The input \code{x} matrix \cr
#'  \code{y} \tab The input \code{y} matrix (\code{\link{Surv}} object)\cr
#'  \code{type} \tab The type of estimation (\code{"lasso"}) \cr
#'  \code{family} \tab The input distribution \cr
#' }
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{relax.lasso}}, \code{ada.lasso}, \code{\link{iregnet}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' 
lasso <- function(x, y, family,...){
    fit <- iregnet(x=x, y=y, family=family,...)
    return(list(coef=fit$beta, scale=fit$scale, lambda=fit$lambda, x=x, y=y, type="lasso", family=family))
}

#' @title Fit relaxed lasso for right-censored AFT model
#' @export
#' @description
#' This is a simple wrapper function to fit accelerated failure time models 
#' using relaxed elastic net penalized maximum likeihood, calling \code{iregnet} and \code{survival::survreg}.
#' Supports gaussian, logistic, Weibull and extreme value distributions.
#'
#' @param x Input matrix of covariates with dimension n_obs * n_vars, with
#' \eqn{nvars \ge 2}. Sparse matrices are not supported.
#' 
#' @param y Response variable. Here, only \code{\link{Surv}} object is supported.
#'
#' @param family The distribution to fit. It can be one of "gaussian", "logistic",
#' "loggaussian", "loglogistic", "extreme_value", "exponential", "weibull". Partial matching
#' is allowed.
#' \cr \emph{Default: "gaussian"}
#'
#' @param ... Further arguments passed to \code{\link{iregnet}} 
#' 
#' 
#' @details This is a wrapper function for \code{\link{iregnet}} with simpler
#' output, facilitating simple cross-validation. The function first calls
#' \code{\link{iregnet}}, and uses the active predictors at each \code{lambda}
#' (the sequence of \code{lambda} is determined by \code{\link{iregnet}})
#' to fit an AFT using \code{\link{survival::survreg}} for each \code{lambda}.
#' @return Returns an object  with the following elements:\cr
#' \tabular{ll}{
#'  \code{coef} \tab Matrix of size \code{(n_vars+1) * num_lambda} containing
#'  intercept, coefficients of \code{X} for each \code{lambda} in the fit model.
#'    \cr
#'  \code{scale} \tab Vector of size \code{num_lambda} of estimated
#'    scale at each \code{lambda} value. \cr
#'  \code{lambda} \tab Vector of size \code{num_lambda} of (calculated or
#'   supplied) regularization parameter \code{lambda} values. \cr
#'  \code{x} \tab The input \code{x} matrix \cr
#'  \code{y} \tab The input \code{y} matrix (\code{\link{Surv}} object)\cr
#'  \code{type} \tab The type of estimation (\code{"relax.lasso"}) \cr
#'  \code{family} \tab The input distribution \cr
#' }
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{lasso}}, \code{ada.lasso}, \code{\link{iregnet}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- relax.lasso(x=X, y=y, family="weibull")
#' 
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


#' @title Fit adaptive lasso for right-censored AFT model
#' @export
#' @description
#' This is a simple wrapper function to fit accelerated failure time models 
#' using adaptive lasso, calling \code{iregnet} and \code{survival::survreg}.
#' Supports gaussian, logistic, Weibull and extreme value distributions.
#'
#' @param x Input matrix of covariates with dimension n_obs * n_vars, with
#' \eqn{nvars \ge 2}. Sparse matrices are not supported.
#' 
#' @param y Response variable. Here, only \code{\link{Surv}} object is supported.
#'
#' @param family The distribution to fit. It can be one of "gaussian", "logistic",
#' "loggaussian", "loglogistic", "extreme_value", "exponential", "weibull". Partial matching
#' is allowed.
#' \cr \emph{Default: "gaussian"}
#' @param omega The $\omega$ parameter of adaptive lasso (default=1)
#' @param ... Further arguments passed to \code{\link{iregnet}} 
#' 
#' 
#' @details This is a wrapper function for \code{\link{iregnet}} with simpler
#' output, facilitating simple cross-validation. The function first calls
#' \code{\link{survival::survreg}}, and then uses the reciprocal absolute regression coefficients
#' raised to the power of $\omega$ as weights of the parameters in a subsequent lasso fit with \code{\link{iregnet}}.
#' As \code{\link{iregnet}} does not support penalty factors (such as \code{glmnet}), this is
#' achieved by multiplying covariates with absolute regression coefficients from the first step
#' and turning off standardization. After the lasso step, the regression coefficients are
#' rescaled. The sequence of \code{lambda} is determined by \code{\link{iregnet}}.
#' @return Returns an object  with the following elements:\cr
#' \tabular{ll}{
#'  \code{coef} \tab Matrix of size \code{(n_vars+1) * num_lambda} containing
#'  intercept, coefficients of \code{X} for each \code{lambda} in the fit model.
#'    \cr
#'  \code{scale} \tab Vector of size \code{num_lambda} of estimated
#'    scale at each \code{lambda} value. \cr
#'  \code{lambda} \tab Vector of size \code{num_lambda} of (calculated or
#'   supplied) regularization parameter \code{lambda} values. \cr
#'  \code{omega} \tab The input \code{omega} \cr
#'  \code{x} \tab The input \code{x} matrix \cr
#'  \code{y} \tab The input \code{y} matrix (\code{\link{Surv}} object)\cr
#'  \code{type} \tab The type of estimation (\code{"ada.lasso"}) \cr
#'  \code{family} \tab The input distribution \cr
#' }
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{relax.lasso}}, \code{lasso}, \code{\link{iregnet}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- ada.lasso(x=X, y=y, family="weibull")
#' 
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

#' @title Plot lasso solution path for right-censored AFT model
#' @export
#' @description
#' This is a simple plot function to visualize the solution path of a lasso,
#' relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{lasso}}, \code{\link{relax.lasso}} or \code{\link{ada.lasso}}.
#' 
#' @param xvar the meaning of the x-axis: either \code{"log"} (default; the logarithm of $\lambda$) 
#' or \code{"L1norm"}, the sum of absolute regression coefficients (except the intercept). 
#' 
#' @details The function plots solution paths against the penalty parameter or the L1 norm of the regularized regression
#' coefficients. It uses standard colors to describe variables 1, 2, 3, ... A dashed line is used for the intercept.
#' @return no return value
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{relax.lasso}}, \code{ada.lasso}, \code{\link{lasso}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' plotlasso(fit)
#' 
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

#' @title Cross-validate lasso, relaxed lasso or adaptive lasso for right-censored AFT model
#' @export
#' @description
#' This is a simple function to crossvalidate a lasso,
#' relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{lasso}}, \code{\link{relax.lasso}} or \code{\link{ada.lasso}}.
#' 
#' @param nfolds the number of folds for cross-validation (default=10)
#' 
#' @param foldid optional: assigned fold id for each observation. If supplied, overrides \code{nfolds}. It is
#' expected that \code{foldid} contains integer numbers of 1 to \code{nfolds}.
#' 
#' @param rep number of repetitions of cross-validation (default=1)
#' 
#' @details The function cross-validates a lasso, relaxed lasso or adaptive lasso to determine the optiomal 
#' value of $\lambda$. Repeated cross-validation is supported.
#' @return Returns an object  with the following elements:\cr
#' \tabular{ll}{
#'  \code{lasso.fit} \tab input fit object \cr
#'  \code{cvm} \tab cross validated log likelihood at each \code{lambda}. The log likelihood is computed using a null fit from 
#'  \code{\link{survival::survreg}} and hence is already normalized to scale of time variable
#'  \code{cvm.sd} \tab standard deviation of fold-specific cross-validated log likelihood at each \code{lambda}.
#'  \code{s} \tab index of optimal \code{lambda}
#'  \code{lambda.min} \tab optimal \code{lambda} (actually the \code{lambda} at which CV log likelihood is maximized)
#'  \code{lambda} \tab Vector of size \code{num_lambda} of (calculated or
#'   supplied) regularization parameter \code{lambda} values. \cr
#'  \code{omega} \tab The input \code{omega} \cr
#'  \code{x} \tab The input \code{x} matrix \cr
#'  \code{y} \tab The input \code{y} matrix (\code{\link{Surv}} object)\cr
#'  \code{type} \tab The type of estimation (taken from \code{fit$type}) \cr
#'  \code{family} \tab The input distribution (taken from \code{fit$family}) \cr
#' }
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{relax.lasso}}, \code{ada.lasso}, \code{\link{lasso}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' 
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

#' @title Plot cross-validated log likelihood vs. lambda for right-censored AFT model
#' @export
#' @description
#' This is a simple plot function to visualize the cross-validation results for a lasso,
#' relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{cv.lasso}}.
#' 
#' @param adj adjustment value to change the spacing between nonzero values of \code{lambda} and 0 (default=5).
#' 
#' @details The function plots cross-validated log likelihood and standard errors against the penalty parameter.
#' The optimal solution is indicated by a horizontal and a vertical line, and a filled dot. The upper x-axis indicates
#' the number of active predictors.
#' @return no return value
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{cv.lasso}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' plot.cv.lasso(cv.fit)
#' 
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

#' @title Optimal coefficients for lasso, relaxed lasso or adaptive lasso fit of right-censored AFT model
#' @export
#' @description
#' This is a simple function to extract the optimal solution as determined by cross-validation from a
#' lasso, relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{cv.lasso}}.
#' 
#' @details The function returns the optimal regression coefficients for lasso, relaxed lasso or 
#' adaptive lasso fit of a right-censored AFT model as determined by cross-validation.
#' @return a named vector of length \code{nvar+1} with the optimal regression coefficients
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{cv.lasso}}, \code{\link{scale.cv.lasso}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' coef.cv.lasso(cv.fit)
#' 
coef.cv.lasso <- function(fit){
    return(fit$lasso.fit$coef[,fit$s])
}

#' @title Optimal coefficients for lasso, relaxed lasso or adaptive lasso fit of right-censored AFT model
#' @export
#' @description
#' This is a simple function to extract the optimal solution as determined by cross-validation from a
#' lasso, relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{cv.lasso}}.
#' 
#' @details The function returns the scale parameter at the optimal regression coefficients for 
#' lasso, relaxed lasso or adaptive lasso fit of a right-censored AFT model as determined by cross-validation.
#' @return numeric: the optimal scale parameter
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{cv.lasso}}, \code{\link{scale.cv.lasso}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' scale.cv.lasso(cv.fit)
#' 
scale.cv.lasso <- function(fit){
    return(fit$lasso.fit$scale[fit$s])
}

#' @title Predictions from lasso, relaxed lasso or adaptive lasso fit of right-censored AFT model
#' @export
#' @description
#' This is a simple function to predict from the optimal solution as determined by cross-validation of a
#' lasso, relaxed lasso or adaptive lasso fit.
#'
#' @param fit object returned by \code{\link{cv.lasso}}.
#' 
#' @param newdata a data set with the active predictors of the optimal solution.
#' 
#' @param type Type of prediction. See documentation of \code{\link{survival::predict.survreg}}. Default=\code{"response"}, which 
#' predicts expected time to failure.
#' 
#' @details The function predicts from the optimal solution as determined by cross-validation of a
#' lasso, relaxed lasso or adaptive lasso fit. It first generates a null fit from \code{survival::survreg} with the optimal solution
#' from \code{cv.lasso}. Hence, it inherits the functionality of \code{survival::predict.survreg}.
#' @return a vector with predicted values for each row of \code{newdata}
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{cv.lasso}}, \code{\link{survival::predict.survreg}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' predict.cv.lasso(cv.fit, newdata=ovarian)
#' 
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

#' @title Survival predictions at specified time from lasso, relaxed lasso or adaptive lasso fit of right-censored AFT model
#' @export
#' @description
#' This is a very simple function to predict survival from the optimal solution as determined by cross-validation of a
#' lasso, relaxed lasso or adaptive lasso fit at a specified time point.
#'
#' @param fit object returned by \code{\link{cv.lasso}}.
#' 
#' @param newdata a data set with the active predictors of the optimal solution.
#' 
#' @param time The time point for which survival should be predicted.
#' 
#' @param grid a grid of survival probabilities
#' 
#' @param se whether standard errors of survival probabilities should be returned (default=\code{TRUE})
#' 
#' @details The function uses a very simple method to determine survival probabilities from the optimal solution as determined by cross-validation of a
#' lasso, relaxed lasso or adaptive lasso fit at a specified time point. As \code{survival::predict.survreg} does not predict
#' survival probabilities but only quantiles, it first determines quantiles of the predicted survival curves for each row of \code{newdata}
#' and then looks up the survival probability which corresponds to a time just below the specified time point. The granularity of \code{grid} 
#' determines the precision of the survival probabilities.
#' @return a vector with predicted survival probabilities for each row of \code{newdata}
#' @author
#' Georg Heinze
#' @useDynLib iregnet
#' @seealso
#' \code{\link{cv.lasso}}, \code{\link{survival::predict.survreg}}
#' @import survival
#' @examples
#' library(survival)
#' X <- cbind(ovarian$ecog.ps, ovarian$rx)
#' y <- Surv(ovarian$futime, ovarian$fustat)
#' fit <- lasso(x=X, y=y, family="weibull")
#' cv.fit <- cv.lasso(fit)
#' surv.cv.lasso(cv.fit, newdata=ovarian, time=)
#' 
surv.cv.lasso <- function(fit, newdata, time, grid=seq(0.001, 0.999, 0.001), se=TRUE){
    quant<-predict.cv.lasso(fit, newdata=newdata, type="quantile", p=grid, se=TRUE)
    return(unlist(sapply(1:nrow(quant$fit), function(X) 1-max(tail(grid[quant$fit[X,]<=time],1),0))))
}

