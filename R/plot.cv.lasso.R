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
#' k<-10
#' n<-300
#' beta <- seq(0,1,1/(k-1))
#' X <- matrix(rnorm(k*n), n, k)
#' failtime <- rexp(n, 1/exp(10 + X %*% beta))
#' maxfu <- quantile(failtime, 0.5)
#' futime <- runif(n, 0, maxfu)
#' status <- (failtime < futime)*1
#' time <- pmin(failtime, futime)
#' y <- Surv(time, status)

#' fit<-lasso(y=y, x=X, family="weibull")
#' cvfit <- cv.lasso(fit)
#' plot.cv.lasso(cvfit)
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
