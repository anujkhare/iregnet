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
#'  \code{call} \tab the call \cr
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
#' @import glmnet
#' @examples
#' library(survival)
#' k<-10
#' n<-300
#' beta <- c(rep(0, 4), seq(0.5, 2, length.out=6))
#' X <- matrix(rnorm(k*n), n, k)
#' failtime <- rexp(n, 1/exp(10 + X %*% beta))
#' maxfu <- quantile(failtime, 0.5)
#' futime <- runif(n, 0, maxfu)
#' status <- (failtime < futime)*1
#' time <- pmin(failtime, futime)
#' y <- Surv(time, status)

#' fit<-lasso(y=y, x=X, family="weibull")
#' plotlasso(fit, xvar="L1norm", intercept=FALSE)
#' 
relax.lasso <- function(x, y, family, ...){
    require(survival)
    if(is.null(colnames(x))) colnames(x)<-paste("X", 1:ncol(x), sep="")
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
    return(list(call=match.call(), coef=coef, scale=scale, lambda=fit$lambda, x=x, y=y, type="relax.lasso", family=family))
}
