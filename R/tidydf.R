#' @title Return a tidy data.frame from iregnet fit
#'
#' @description
#' Returns a tidy \link{data.frame} from a fitted "iregnet" object.
#'
#' @param x The S3 object of type \code{iregnet} returned by the \code{iregnet}
#' method.
#'
#' @param ... Other parameters. Currently unused.
#'
#' @details
#' This function is used to obtain an intermediate \code{data.frame} used in
#' \link{plot.iregnet}.
#' It can be used for producing other plots using \code{ggplot2}.
#' NOTE: \code{Intercept} (if present) is \strong{not} included in the
#' \code{arclength} since it is never regularized.
#'
tidydf <- function(x, ...) {
  stopifnot_error("Invalid / no fit object provided",
                  class(x) == "iregnet")

  # Don't include intercept in norm (arclength) since it is never regularized.
  start_index <- as.integer(x$intercept) + 1
  n <- nrow(x$beta)
  arclength <- apply(x$beta, 2, function(x) sum(abs(x[start_index: n])))
  tidy.df <- with(x, data.frame(
    weight=as.numeric(t(beta)),
    lambda,
    arclength,
    variable=rownames(beta)[as.integer(col(t(beta)))]
    ))
  tidy.df
}
