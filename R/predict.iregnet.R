#' @title Predict response using new covariates
#'
#' @description
#' Prediction is just X'beta for non-transformed distributions, and
#' itrans(X'beta) for transformed distributions.
#' In the case of \code{gaussian}, \code{logistic} and
#' \code{extreme_value}, types link and response both return the linear
#' predictors.
#'
#' @param fit The S3 object of type \code{iregnet} returned by the \link{iregnet}
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
#' \cr \emph{Default: \code{NULL} (\code{fit$lambda} is used)}
#'
#' @param ... Optional arguments. Currently unused.
predict.iregnet <- function(fit, newx, lambda=NULL, type=c("link", "response"), ...) {
  stopifnot_error("Invalid / no fit object provided", !missing(fit),
    class(fit) == "iregnet")
  stopifnot_error("No 'newx' matrix provided", !missing(newx))
  stopifnot_error("newx should be a matrix with same number of columns as used in fit",
    is.matrix(newx), ncol(newx) + fit$intercept == nrow(fit$beta))

  # ensure that all the lambda values are in the fit
  if (is.null(lambda))
    lambda = fit$lambda
  inds <- match(lambda, fit$lambda)
  stopifnot_error("Lambda values must be those used in fit", all(!is.na(inds)))

  # append intercept column
  if (fit$intercept)
    newx <- cbind2(1, newx)

  beta <- fit$beta[, inds]
  link <- newx %*% beta  # link == eta
  response <- link

  if (fit$family %in% names(transformed_distributions)) {
    response <- transformed_distributions[[fit$family]]$trans(link)
  }

  type <- match.arg(type)
  switch(type,
         link = link,
         response = response)
}  
