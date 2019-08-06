#' GEV generalized linear regression modelling
#'
#' Regression modelling with a Generalized Extreme value response distribution
#' and generalized linear modelling of each parameter, specified by a symbolic
#' description of each linear predictor and the inverse link function to be
#' applied to each linear predictor.  Models are fitted using maximum
#' likelihood estimation using the \code{\link[stats]{optim}} function.
#' @param y Either a numeric vector or the name of a variable in \code{data}.
#'   \code{y} must not have any missing values.
#' @param data An optional data frame containing \code{y} and any covariates.
#'   If these variables are not found in \code{data} then the variables are
#'   taken from the environment from which \code{gevreg} is called.
#'   Neither \code{y} nor the covariates may have any missing values.
#' @param mu,sigma,xi Formulae (see \code{\link[stats]{formula}})
#'   for the GEV parameters \code{mu} (location), \code{sigma} (scale)
#'   and \code{xi} (shape), e.g. \code{mu = ~ x}.
#' @param mustart,sigmastart,xistart Optional numeric vectors specifying
#'   respective initial values for the parameters relating to location,
#'   scale and shape.  If not supplied these are set inside the
#'   \code{gevreg}.
#' @param mulink,sigmalink,xilink Functions giving the respective
#'   link functions that relate the location, scale and
#'   shape parameters to the linear predictor, for example
#'   \code{sigmalink(sigma) =} linear predictor.
#'   The inverses of these functions may be supplied using
#'   \code{invmulink, invsigmalink, invxilink}, but, otherwise, if a link
#'   is one of \code{identity} \code{log} then the
#'   corresponding inverse link function is from the link function.
#'   Otherwise, these functions are only
#'   used to set initial estimates of the parameters. \code{mulink} is not
#'   used if \code{mustart} is supplied, etc.
#' @param invmulink,invsigmalink,invxilink Functions giving the respective
#'   \strong{inverse} link functions that relate the location, scale and
#'   shape parameters to the linear predictor.  If these are supplied then
#'   the code does \strong{not} check that they are consistent with the link
#'   functions supplied in \code{mulink,sigmalink,xilink}.
#' @param optim_control A list to be passed to \code{\link[stats]{optim}}
#'   as its argument \code{control}.
#' @param ... further arguments to be passed to \code{\link[stats]{optim}}.
#' @details Add details.
#' @warning If the input data contains covariate 'year', 'year' has to be
#'   scaled before fitting it.
#' @return An object (list) of class \code{c("gev", "evreg")}, which has
#'   the following components
#'     \item{coefficients}{A named numeric vector of the estimates of the
#'       model parameters.}
#'     \item{se}{estimated standard errors}
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @examples
#' ### Oxford-Worthing temperature data
#'
#' ## Estimate separate GEV parameters for for Oxford and Worthing,
#' ## as in Chandler and Bate (2007)
#' # Intercepts are Oxford-Worthing average
#' ow1 <- gevreg(temp, data = ow, mu = ~ loc, sigma = ~ loc,
#'               xi = ~ loc, sigmalink = identity)
#' # Intercepts relate to Oxford
#' ow2 <- gevreg(temp, data = ow, mu = ~ factor(name), sigma = ~ factor(name),
#'               xi = ~ factor(name), sigmalink = identity)
#'
#' ### Fremantle sea levels
#'
#' ## No covariates
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' ## Add SOI as a covariate
#' f1 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ SOI)
#' ## Add (shited and scaled) year as a covariate (instead)
#' f2 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year01)
#' # Include both SOI and year
#' f3 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year01 + SOI)
#' # Example of user-supplied link
#' f3 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year01 + SOI,
#'              mulink = function(x) x, invmulink = function(x) x)
#' @export
gevreg <- function(y, data, mu = ~ 1, sigma = ~ 1, xi = ~ 1,
                   mustart, sigmastart, xistart, mulink = identity,
                   sigmalink = log, xilink = identity, invmulink,
                   invsigmalink, invxilink,
                   optim_control = list(maxit = 10000), ...) {
  # Record the call for later use
  Call <- match.call()
  # Check that the data have been supplied
  if (missing(data)) {
    stop("data must be supplied")
    data <- NULL
  } else {
    y <- deparse(substitute(y))
  }
  # Create intercept formula for every parameter
  gev_pars <- c("mu", "sigma", "xi")
  mp <- c(mu = mu, sigma = sigma, xi = xi)
  # Create the data
  model_data <- create_data(y = y, data = data, params = mp)
  # Infer missing inverse link functions from link function
  if (missing(invmulink)) {
    namelink <- deparse(substitute(mulink))
    invmulink <- infer_invlinks(namelink)
  }
  if (missing(invsigmalink)) {
    namelink <- deparse(substitute(sigmalink))
    invsigmalink <- infer_invlinks(namelink)
  }
  if (missing(invxilink)) {
    namelink <- deparse(substitute(xilink))
    invxilink <- infer_invlinks(namelink)
  }
  # Set starting values, if necessary
  if (missing(mustart)) {
    mustart <- gevreg_mustart(model_data, mulink)
  }
  if (missing(sigmastart)) {
    sigmastart <- gevreg_sigmastart(model_data, sigmalink)
  }
  if (missing(xistart)) {
    xistart <- gevreg_xistart(model_data, xilink)
  }
  start <- c(mustart, sigmastart, xistart)

  # Estimate the parameters using stats::optim
  res <- stats::optim(par = start, fn = gevreg_negloglik, data = model_data,
                      invmulink = invmulink, invsigmalink = invsigmalink,
                      invxilink = invxilink, hessian = TRUE,
                      control = optim_control, ...)
  # Add information to the returned object
  res$call <- Call
  res$data <- model_data
  res$coefficients <- res$par
  nms <- unlist(lapply(names(res$data$D),
                       function(x){
                         paste(x, ": ", colnames(res$data$D[[x]]), sep = "")
                         }))
  names(res$coefficients) <- nms
  res$formulae <- mp
  res$links <- list(mulink = mulink, sigmalink = sigmalink, xilink = xilink)
  res$invlinks <- list(invmulink = invmulink, invsigmalink = invsigmalink,
                       invxilink = invxilink)
  res$residuals <- gevreg_residuals(res)
  res$loglik <- -res$value
  res$cov <- solve(res$hessian)
  res$se <- sqrt(diag(res$cov))
  res$value <- res$counts <- res$hessian <- NULL
  res$df <- length(res$coefficients)
  # We don't need xlevels, for the moment
#  res$xlevels <- gevreg_xlevels(res$formulae, data)
  class(res) <- c("gev", "evreg")
  res$aic <- AIC(res)
  return(res)
}
