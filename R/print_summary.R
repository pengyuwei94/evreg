# ================================ print.evreg ================================

#' Print method for objects of class "evreg"
#'
#' \code{print} method for class "evreg".
#'
#' @param x an object inheriting from class "evreg", a result of a call to
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
    stop("use only with \"evreg\" objects")
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

# ============================== summary.evreg ================================

#' Summarizing GEV regression model fits
#'
#' These functions are all methods for class \code{evreg} or
#' \code{summary.evreg} objects.  They provide similar functionality
#' to \code{\link[stats]{summary.glm}}.
#'
#' @param object an object inheriting from class "evreg", a result of a call to
#'   \code{\link{gevreg}}.
#' @param correlation logical; if TRUE, the correlation matrix of the estimated
#'   parameters is returned and printed.
#' @param symbolic.cor logical. If TRUE, print the correlations in a symbolic
#'   form (see \code{\link[stats]{symnum}}) rather than as numbers.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{format}} and \code{\link[base:Round]{signif}}.
#' @param signif.stars logical. If TRUE, 'significance stars' are printed for
#'   each coefficient.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return \code{summary.glm} returns an object of class
#'   \code{"summary.evreg"}, a list with components
#'   \item{call}{the component \code{call} from \code{object}.}
#'   \item{coefficients}{the matrix of coefficients, standard errors, z-values
#'     and p-values.}
#'   \item{deviance}{the deviance of the fitted model, calculated as
#'     \code{-2 * logLik(object)}}
#'   \item{df.residual}{the residual degrees of freedom, calculated as
#'     \code{length(object$data$y) - object$df}}
#'   \item{aic}{the component \code{aic} from \code{object}.}
#'
#'   \code{print.summary.glm} prints the call, a table of estimated
#'   coefficients, standard errors, the estimate/SE ratio and a two-tailed
#'   approximate p-value based on a standad normal reference distribution.
#'   Significance stars are added if \code{signif.stars = TRUE}.
#'   Also printed are the residual deviance and residual degrees of freedom
#'   and the value of AIC.
#'   The argument \code{x} is returned invisibly.

#' @seealso \code{\link{gevreg}} GEV generalized linear regression modelling.
#' @name summary.evreg
NULL
## NULL

#' @rdname summary.evreg
#' @export
summary.evreg <- function(object, correlation = FALSE, symbolic.cor = FALSE,
                          ...) {
  if (!inherits(object, "evreg")) {
    stop("use only with \"evreg\" objects")
  }
  res <- list()
  res$call <- object$call
  coefs <- coef(object)
  zvalue <- coef(object) / object$se
  pvalue <- 2 * stats::pnorm(-abs(zvalue))
  coef_table <- cbind(coefs, object$se, zvalue, pvalue)
  dimnames(coef_table) <- list(names(coefs), c("Estimate", "Std. Error",
                                               "z value", "Pr(>|z|)"))
  res$coefficients <- coef_table
  if (correlation) {
    res$correlation <- stats::cov2cor(vcov(object))
    res$symbolic.cor <- symbolic.cor
  }
  res$deviance <- -2 * logLik(object)
  res$df.residual <- length(object$data$y) - object$df
  res$aic <- object$aic
  res$pvalue <- pvalue
  class(res) <- "summary.evreg"
  return(res)
}

# ============================= print.summary.evreg ===========================

#' @rdname summary.evreg
#' @export
print.summary.evreg <- function (x, digits = max(5L, getOption("digits") - 3L),
                                 symbolic.cor = x$symbolic.cor,
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  if (!inherits(x, "summary.evreg")) {
    stop("use only with \"summary.evreg\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, na.print = "NA", ...)
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- ncol(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\nResidual deviance:", format(signif(x$deviance, digits)),
      " on", format(x$df.residual), " degrees of freedom")
  cat("\nAIC:", format(signif(x$aic, digits)))
  invisible(x)
}
