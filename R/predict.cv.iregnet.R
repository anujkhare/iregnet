##' Predict CV iregnet
##'
##' Get predicted values
##' @title predict CV iregnet
##' @param object result of cv.iregnet
##' @param newx feature matrix
##' @param type min or 1sd
##' @param ... further parameters
##' @return predicted values
##' @export
##' @author Toby Dylan Hocking
predict.cv.iregnet <- function(object, newx, type="min", ...){
  i <- object$selected[[type]]
  not.cv <- object
  class(not.cv) <- "iregnet"
  predict(not.cv, newx, lambda=object$lambda[[i]], ...)
}