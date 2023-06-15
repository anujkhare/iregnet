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
