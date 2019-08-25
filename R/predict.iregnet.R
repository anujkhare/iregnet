#' @title Predict response using new covariates
#' @export
#' @description
#' Prediction is just X'beta for non-transformed distributions, and
#' itrans(X'beta) for transformed distributions.
#' In the case of \code{gaussian}, \code{logistic} and
#' \code{extreme_value}, types link and response both return the linear
#' predictors.
#'
#' @param object The S3 object of type \code{iregnet} returned by the \link{iregnet}
#' method.
#'
#' @param newx New input matrix of covariates with dimension n_vars same as used
#' for fit. Sparse matrices are not supported.
#'
#' @param type Type of prediction required. Type \code{link} returns the linear
#' predictors. Type \code{response} returns the fitted values.
#' In the case of \code{gaussian}, \code{logistic} and
#' \code{extreme_value}, the link and response are the same.
#'
#' @param lambda The values of lambda at which prediction is to be obtained.
#' These should be a subset of the values on which the model was fit. To obtain
#' predictions at other lambda values, re-fit the model.
#' \cr \emph{Default: \code{NULL} (\code{object$lambda} is used)}
#' @import methods
#' @import stats
#' @param ... Optional arguments. Currently unused.
predict.iregnet <- function(object, newx, lambda=NULL, type=c("link", "response"), ...) {
  stopifnot_error("No 'newx' matrix provided", !missing(newx))
  # ensure that all the lambda values are in the fit
  if (is.null(lambda))
    lambda = object$lambda
  inds <- match(lambda, object$lambda)
  stopifnot_error("Lambda values must be those used in fit", all(!is.na(inds)))
  # check non-zero 
  beta <- object$beta[, inds, drop=FALSE]
  is.used <- apply(beta!=0, 1, any)
  not.intercept <- rownames(beta) != "(Intercept)"
  feature.name.vec <- rownames(beta)[is.used & not.intercept]
  stopifnot_error("newx should be a numeric matrix",
                  is.matrix(newx), is.numeric(newx))
  has.feature <- feature.name.vec %in% colnames(newx)
  feature.not.present <- feature.name.vec[!has.feature]
  if(length(feature.not.present))
    colnames(newx) <- paste('x', 1: ncol(newx), sep='')
  # if(length(feature.not.present) == nrow(newx))
  #   colnames(newx) <- paste('x', 1: ncol(newx), sep='')
  # else if(length(feature.not.present)){
  #   stop("features missing but needed for prediction: ",
  #        paste(feature.not.present, collapse=", "))
  # }

  one <- if(object$intercept)1
  needed.features <- cbind(one, newx[, feature.name.vec, drop=FALSE])
  coef.name.vec <- c(if(object$intercept)"(Intercept)", feature.name.vec)
  needed.beta <- beta[coef.name.vec, , drop=FALSE]
  link <- needed.features %*% needed.beta  # link == eta
  response <- if (object$family %in% names(transformed_distributions)) {
    transformed_distributions[[object$family]]$trans(link)
  }else{
    link
  }
  type <- match.arg(type)
  switch(type,
         link = link,
         response = response)
}  
