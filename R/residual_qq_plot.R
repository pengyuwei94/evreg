# ================================= plot.evreg ============================== #

#' Plot method for evreg fitted model objects
#'
#' Produces a residual QQ plot based on fitted model objects returned from
#' \code{\link{gevreg}} and \code{\link{ppreg}}.
#'
#' @param x An object returned from either \code{\link{gevreg}} or
#'   \code{\link{ppreg}}.
#' @param y Not used.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}}.
#' @details If \code{x} was returned from \code{\link{gevreg}} then the
#'  ordered (smallest to largest) residuals \eqn{r_i, i = 1, ..., n} are
#'  plotted against the 100\eqn{i / (n + 1)}\% \eqn{i = 1, ..., n}
#'  quantiles of a standard Gumbel distribution.
#'  If \code{x} was returned from \code{\link{ppreg}} then the
#'  ordered residuals are plotted against the
#'  100\eqn{i / (n + 1)}\% \eqn{i = 1, ..., n} of a standard exponential
#'  distribution.
#' @return Nothing is returned, only the plot is produced.
#' @examples
#' # GEV
#' f1 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ SOI)
#' plot(f1)
#'
#' # PP (no coavriates)
#' pp <- ppreg(y = SeaLevel, data = fremantle[, -1])
#' plot(pp)
#'
#' # PP regression
#' pp <- ppreg(y = SeaLevel, data = fremantle[, -1], mu = ~SOI)
#' @export
plot.evreg <- function(x, y = NULL, ...) {
  if (!inherits(x, "evreg")) {
    stop("use only with \"evreg\" objects")
  }
  model <- class(x)[1]
  n <- length(x$residuals)
  xvals <- (1:n) / (n + 1)
  my_ylab <- "ordered residuals"
  if (model == "gev") {
    my_xlab <- "Gumbel quantiles"
    xvals <- -log(-log(xvals))
  } else if (model == "pp") {
    my_xlab <- "exponential quantiles"
    xvals <- -log(1 - xvals)
  }
  my_plot <- function(x, y, ..., xlab = my_xlab, ylab = my_ylab) {
    graphics::plot(x = x, y = y, ..., xlab = xlab, ylab = ylab)
  }
  my_plot(x = xvals, y = sort(x$residuals), ...)
  graphics::abline(0, 1, lty = 1, col = "blue")
  if (model == "gev") {
    my_title <- "Residual QQ Plot (Gumbel Scale)"
  } else if (model == "pp") {
    my_title <- "Residual QQ Plot (Exponential Scale)"
  }
  graphics::title(main = my_title)
  invisible()
}
