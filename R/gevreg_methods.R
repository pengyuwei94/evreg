#' @export
coef.evreg <- function(object, ...){
  return(object$coefficients)
}

#' @export
nobs.evreg <- function(object, ...) {
  return(length(object$data$y))
}

#' @export
vcov.evreg <- function(object, ...) {
  vc <- object$cov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.evreg <- function(object, ...) {
  val <- object$loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
