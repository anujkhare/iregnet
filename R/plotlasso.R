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
#' @param intercept if TRUE (default), adds the intercept path to the plot
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
plotlasso <- function(fit, xvar="log", intercept=TRUE){
    if(xvar=="log") {
        xv <- log(fit$lambda+1)
        xlabel <- "log(lambda+1)"
    } else if(xvar=="L1norm") {
        xv <- apply(fit$coef[-1,],2,function(X) sum(abs(X)))
        xlabel <- "L1 norm"
    }
    if(!intercept) coef <- fit$coef[-1,]
    else coef <- fit$coef
    ylim <- range(coef)
    xlim <- range(xv)
    plot(x=xv, y=fit$coef[1,], type="l", lty=3, xlim=xlim, ylim=ylim, ylab="beta", xlab=xlabel)
    for(j in 1:nrow(fit$coef[-1,])) lines(x=xv, y=fit$coef[j+1,], type="l", lty=1, col=j)
}
