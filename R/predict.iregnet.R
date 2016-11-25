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
  stopifnot_error("Invalid / no fit object provided", !missing(object),
    class(object) == "iregnet")
  stopifnot_error("No 'newx' matrix provided", !missing(newx))
  stopifnot_error("newx should be a matrix with same number of columns as used in fit",
    is.matrix(newx), ncol(newx) + object$intercept == nrow(object$beta))

  # ensure that all the lambda values are in the fit
  if (is.null(lambda))
    lambda = object$lambda
  inds <- match(lambda, object$lambda)
  stopifnot_error("Lambda values must be those used in fit", all(!is.na(inds)))

  # append intercept column
  if (object$intercept)
    newx <- cbind2(1, newx)

  beta <- object$beta[, inds]
  link <- newx %*% beta  # link == eta
  response <- link

  if (object$family %in% names(transformed_distributions)) {
    response <- transformed_distributions[[object$family]]$trans(link)
  }

  type <- match.arg(type)
  switch(type,
         link = link,
         response = response)
}  
