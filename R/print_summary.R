# ================================ print.evreg ================================

#' Print method for objects of class "evreg"
#'
#' \code{print} method for class "evreg".
#'
#' @param x an object of class "evreg", a result of a call to
#'   \code{\link{gevreg}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{format}} and \code{\link[base:Round]{signif}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details Prints: the call, an indicator of convergence (\code{TRUE} only
#'   if \code{x$call = 0}), estimated coefficients, maximized log-likelihood,
#'   the number of estimated coefficients (DF) and AIC.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{gevreg}} GEV generalized linear regression modelling.
#' @export
print.evreg <- function(x, digits = max(3, getOption("digits") - 3L), ...) {
  if (!inherits(x, "evreg")) {
    stop("use only with \"chandwich\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (x$conv == 0) {
    conv <- TRUE
  } else {
    conv <- FALSE
  }
  cat( "Convergence:", conv, "\n\n")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nLog-likelihood:", format(signif(logLik(x), digits)),
      "\tDF:", x$df,
      "\tAIC:", format(signif(x$aic, digits)))
  return(invisible(x))
}
