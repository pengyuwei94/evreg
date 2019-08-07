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
#'   is one of \code{identity} or \code{log} then the
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
#'   as its argument \code{control}.  The default setting in \code{gevreg}
#'   sets \code{maxit} to 10000 because the default in
#'   \code{\link[stats]{optim}} (100) may not be sufficient when there are
#'   several covariates and/or the initial estimates are poor.
#' @param scale_covs A logical scalar.  Should we center and scale the
#'   covariate data before minmizing the negated log-likelihood using
#'   \code{\link[stats]{optim}}?  Doing this is advisable if the covariate
#'   data have very different orders of magnitude.  The results (parameter
#'   estimates and their estimated covariance matrix) are reported on the
#'   original scale either way.  The fits using \code{scale_covs = TRUE}
#'   and \code{scale_covs = FALSE} may be slightly different.
#'   \code{gevreg} sets \code{optim_control$reltol = 1e-16} to force the fits
#'   to be closer than they would be otherwise.
#' @param ... further arguments to be passed to \code{\link[stats]{optim}}.
#' @details Add details.
#'
#' The default for \code{scale_covs} is \code{FALSE}.  A classic case when a
#' model fit may fail under this setting is when calendar year is a covariate.
#' This covariate will tend to be orders of magnitude larger than other
#' covariates, which causes numerical optimisation problems.  A solution is
#' to use \code{scale_covs = TRUE}.  However, it may be better to shift/scale
#' the calendar year prior to calling \code{gevreg} in order to give the
#' intercept of the model a more useful interpretation.  Otherwise, the
#' intercept corresponds to the year 0, which is unlikely to be of interest.
#'
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
#' summary(f3)
#' f3b <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year01 + SOI,
#'               scale_covs = TRUE)
#' summary(f3b)
#' f3c <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year + SOI,
#'               scale_covs = TRUE)
#' summary(f3c)
#' # (Note: this fit fails if scale_covs = FALSE)
#'
#' # Example of user-supplied link
#' f3 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ Year01 + SOI,
#'              mulink = function(x) x, invmulink = function(x) x)
#' @export
gevreg <- function(y, data, mu = ~ 1, sigma = ~ 1, xi = ~ 1,
                   mustart, sigmastart, xistart, mulink = identity,
                   sigmalink = log, xilink = identity, invmulink,
                   invsigmalink, invxilink,
                   optim_control = list(maxit = 10000, reltol = 1e-16),
                   scale_covs = FALSE, ...) {
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
  if (scale_covs) {
    temp <- scale_covariates(y = y, data = data, mp = mp)
    model_data <- temp$model_data
    Tmat <- temp$Tmat
  } else {
    model_data <- create_data(y = y, data = data, params = mp)
  }
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
  # If necessary, reverse the effects of scaling the covariates
  if (scale_covs) {
    saved_names <- names(res$coefficients)
    res$coefficients <- as.vector(res$coefficients %*% t(Tmat))
    names(res$coefficients) <- saved_names
    res$cov <- Tmat %*% res$cov %*% t(Tmat)
    dimnames(res$cov) <- list(NULL, NULL)
  }
  #
  res$se <- sqrt(diag(res$cov))
  res$value <- res$counts <- res$hessian <- NULL
  res$df <- length(res$coefficients)
  # We don't need xlevels, for the moment
#  res$xlevels <- gevreg_xlevels(res$formulae, data)
  class(res) <- c("gev", "evreg")
  res$aic <- AIC(res)
  return(res)
}

find_Tmat <- function(dmat) {
  mu_data <- dmat
  which_cov <- which(colnames(mu_data) != "(Intercept)")
  # Scale each covariate to have mean 0 and SD 1
  temp <- scale(mu_data[, which_cov])
  dmat[, which_cov] <- temp
  # Extract the centering and scaling values
  cvals <- attributes(temp)$`scaled:center`
  svals <- attributes(temp)$`scaled:scale`
  Tmat1 <- c("(Intercept)" = 1, -cvals / svals)
  diag_mat <- diag(1 / svals)
  zero_mat <- matrix(0, nrow = length(svals), ncol = 1)
  Tmat <- rbind(Tmat1, cbind(zero_mat, diag_mat))
  return(list(dmat = dmat, Tmat = Tmat))
}

scale_covariates <- function(y, data, mp) {
  # When I create Tmat below the covariates must be in the same order as they
  # are in model_data$D$mu, model_data$D$sigma and model_data$D$xi.
  # Call create_data() first to find out which covariates are present in mu,
  # sigma and xi and to get the automatically in the correct order.
  model_data <- create_data(y = y, data = data, params = mp)
  temp <- find_Tmat(model_data$D$mu)
  model_data$D$mu <- temp$dmat
  Tmu <- temp$Tmat
  temp <- find_Tmat(model_data$D$sigma)
  model_data$D$sigma <- temp$dmat
  Tsigma <- temp$Tmat
  temp <- find_Tmat(model_data$D$xi)
  model_data$D$xi <- temp$dmat
  Txi <- temp$Tmat
  # Create the matrix T
  n_mu <- ncol(model_data$D$mu)
  n_sigma <- ncol(model_data$D$sigma)
  n_xi <- ncol(model_data$D$xi)
  zero_sigma_xi <- matrix(0, nrow = n_sigma, ncol = n_xi)
  Tsigma_xi <- rbind(cbind(Tsigma, zero_sigma_xi),
                     cbind(t(zero_sigma_xi), Txi))
  zero_mu <- matrix(0, n_mu, n_sigma + n_xi)
  Tmat <- rbind(cbind(Tmu, zero_mu), cbind(t(zero_mu), Tsigma_xi))
  # Return the transformed data and the parameter transformation matrix
  return(list(model_data = model_data, Tmat = Tmat))
}
