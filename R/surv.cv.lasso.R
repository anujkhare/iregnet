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
