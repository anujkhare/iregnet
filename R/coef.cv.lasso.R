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
