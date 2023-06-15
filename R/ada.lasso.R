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
#'  \code{call} \tab the call \cr
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

#' fit<-ada.lasso(y=y, x=X, family="weibull")
#' plotlasso(fit, xvar="L1norm", intercept=FALSE)
#' 
ada.lasso <- function(x, y, family, omega=1, ...){
    require(survival)
    if(is.null(colnames(x))) colnames(x)<-paste("X", 1:ncol(x), sep="")
    fit_init <- survreg(y~x, dist=family)
    weights <- abs(coef(fit_init)[-1])**omega
    xw <- x %*% diag(weights)
    fit_lasso <- iregnet(y=y, x=xw, standardize=FALSE, family=family, ...)
    ada.beta<-fit_lasso$beta
    ada.beta <- rbind(ada.beta[1,], diag(weights) %*% fit_lasso$beta[-1,])
    rownames(ada.beta) <- c("(Intercept)",colnames(x))
    return(list(call=match.call(), coef=ada.beta, scale=fit_lasso$scale, lambda=fit_lasso$lambda, omega=omega, x=x, y=y, type="ada.lasso", family=family))
}
