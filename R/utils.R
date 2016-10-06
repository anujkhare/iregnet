#' @title stopifnot with custom error message
#'
#' @description
#' Like \code{stopifnot}, but with a custom error message.
#'
#' @param err_message The error message to print.
#' @param ... An error is raised if any these expressions is \code{FALSE}.
stopifnot_error <- function(err_message, ...)
{
  n <- length(ll <- list(...))
  for (i in 1:n)
    if (!(is.logical(r <- ll[[i]]) && !anyNA(r) && all(r))) {
      stop(err_message)
    }
}

#' @title Get censoring status from a \code{Surv} object
#'
#' @description
#' Returns the status column from a \link{Surv} object, and converts to the
#' form: right censored - 0, none - 1, left - 2, interval - 3.
#'
#' @param s \code{Surv} object.
get_status_from_surv <- function(s)
{
  type <- attr(s, 'type')

  stopifnot_error("Unsupported censoring type from Surv", type %in% c('left', 'right',
                                                                      'interval', 'interval2'))
  # right censored: 0, none: 1, left: 2, interval: 3
  status <- s[, ncol(s)]
  if (type == 'left')
      status[status == 0] <- 2

  return(status)
}

#' @title List of non-basic distributions in terms of basic distributions
#'
#' @description
#' A \code{list} containing the supported distributions that are
#' transformations of other "base" distributions. For instance,
#' \code{loggaussian} is the same as \code{gaussian} after log transforming the
#' output \code{y}.
#' @details
#' The list contains the transformation \code{trans}, the base distribution
#' \code{dist},  and the inverse transformation \code{itrans}. It is used
#' internally to fit such distributions.
transformed_distributions <- list(
  "loggaussian" = list(trans = function(y) log(y), itrans = function(y) exp(y), dist = 'gaussian'),
  "loglogistic" = list(trans = function(y) log(y), itrans = function(y) exp(y), dist = 'logistic'),
  "weibull" = list(trans = function(y) log(y), itrans = function(y) exp(y), dist = 'extreme_value'),
  "exponential" = list(trans = function(y) log(y), itrans = function(y) exp(y), dist = 'extreme_value')
)
