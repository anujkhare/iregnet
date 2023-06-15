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
                        data=x,dist=fit$family, init=coef.cv.lasso(fit)[c(TRUE,active)], scale=scale.cv.lasso(fit), control=survreg.control(maxiter=0))
    }
    p<-predict(fit0, newdata, type, ...)
    return(p)
}
