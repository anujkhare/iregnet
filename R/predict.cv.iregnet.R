#' Predict CV iregnet
#'
#' Get predicted values
#' @title predict CV iregnet
#' @param object result of cv.iregnet
#' @param newx feature matrix
#' @param type Type of prediction required. Type \code{link} returns the linear
#' predictors. Type \code{response} returns the fitted values.
#' In the case of \code{gaussian}, \code{logistic} and
#' \code{extreme_value}, the link and response are the same.
#' @param lambda.type min or 1sd
#' @param ... further parameters
#' @return predicted values
#' @export
#' @author Toby Dylan Hocking
predict.cv.iregnet <- function(object, newx, type=c("link", "response"), lambda.type="min", ...){
  i <- object$selected[[lambda.type]]
  not.cv <- object
  class(not.cv) <- "iregnet"
  predict(not.cv, newx, type, lambda=object$lambda[[i]], ...)
}