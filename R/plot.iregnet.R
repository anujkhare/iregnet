#' @title Plot the path of the variables of iregnet fit
#'
#' @description
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' "iregnet" object.
#'
#' @param fit The S3 object of type \code{iregnet} returned by the \code{iregnet}
#' method.
#'
#' @param xvar Variable on the X-axis against which the coefficients are
#' plotted. \code{norm} plots against the L1 norm of the coefficients, i.e.,
#' arclength. \code{lambda} plots against log-lambda sequence.
#' \cr \emph{Default: "norm"}
#'
#' @param label If \code{TRUE}, coefficient names / variable sequence numbers
#' are plotted along with the curves.
#'
#' @param ... Other parameters. Currently unused.
#'
#' @details
#' This function uses \link{tidydf} function to obtain a \code{data.frame}
#' from the \link{iregnet} object. \link{ggplot} is used for the plotting
#' It can be directly used to produce other plots.
#' \code{Intercept} (if present) is \strong{not} included in the \code{arclength}
#' since it is never regularized. It is also not plotted.
#'
plot.iregnet <- function(fit, xvar=c("norm", "lambda"), label=T, ...) {
  stopifnot_error("Invalid / no fit object provided", !missing(fit),
                  class(fit) == "iregnet")
  xvar <- match.arg(xvar)

  tidy.df <- tidydf(fit)
  start_index <- as.integer(fit$intercept) + 1
  n <- nrow(fit$beta)
  varnames <- rownames(fit$beta)[start_index:n]
  tidy.df <- tidy.df[tidy.df$variable %in% varnames, ]
  switch(xvar,
    "lambda" = {
      fig.iregnet.profile <- ggplot()+
        geom_line(aes(log(lambda), weight, color=variable),
                   data=tidy.df)+
        xlab("Log Lambda")+
        ylab("Coefficients")
    },
    "norm" = {
      fig.iregnet.profile <- ggplot()+
        geom_line(aes(arclength, weight, color=variable),
                   data=tidy.df)+
        xlab("L1 Norm of Coefficients")+
        ylab("Coefficients")
    }
  )
  print(fig.iregnet.profile)
}
