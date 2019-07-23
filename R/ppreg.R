#' Point process generalized linear regression modelling
#'
#' Regression modelling with a Point Process response distribution
#' and generalized linear modelling of each parameter, specified by a symbolic
#' description of each linear predictor and the inverse link function to be
#' applied to each linear predictor.
#' Quantile regression has been used to set a threshold for which the
#' probability \code{p} of threshold exceedance is approximatly constant
#' across different values of covarites.
#'
#' @param y Either a numeric vector or the name of a variable in \code{data}.
#'   \code{y} must not have any missing values.
#' @param data Al data frame containing \code{y} and any covariates.
#'   Neither \code{y} nor the covariates may have any missing values.
#' @param p Probability quantile. This is generally a number strictly between 0 and 1.
#' @param mu,sigma,xi Formulae (see \code{\link[stats]{formula}})
#'   for the PP parameters \code{mu} (location), \code{sigma} (scale),
#'   e.g. \code{mu = ~ x}. Parameter \code{xi} (shape) is holded as fixed.
#'
#' @Warning At the moment, only models with identity inverse link function for
#' all parameters and constant shape parameter can be fitted.
#' Different optimization methods may result in wildly different parameter estimates.





ppreg <- function(y, data, p = 0.5, npy = 365, mu = ~1, sigma = ~1, xi = ~1,
                  mustart, sigmastart, xistart, invmulink = identity,
                  invsigmalink = identity, invxilink = identity, ...){
  # Record the call for later use
  Call <- match.call()
  # 1. Check that the data have been supplied
  if (missing(data)) {
    stop("data must be supplied")
    data <- NULL
  } else {
    y <- deparse(substitute(y))
  }

  # 2. Check if the input inverse link functions are identity
  if ((invmulink != identity) ||
      (invsigmalink != identity) || (invxilink != identity)) {
    stop("Inverse link function has to be identity for each parameter")
  }

  # 3. Check if input p is valid
  if(p > 1 || p < 0){
    stop("Probabiliry quantile should be smaller than 1 and larger than 0")
  }


  # Create intercept formula for every parameter
  gev_pars <- c("mu", "sigma", "xi")
  mp <- c(mu = mu, sigma = sigma, xi = xi)
  # Create the data
  model_data <- create_data(y = y, data = data, params = mp)

  # Set starting values, if necessary
  if (missing(mustart)) {
    mustart <- gevreg_mustart(model_data, invmulink)
  }
  if (missing(sigmastart)) {
    sigmastart <- gevreg_sigmastart(model_data, invsigmalink)
  }
  if (missing(xistart)) {
    xistart <- gevreg_xistart(model_data, invxilink)
  }
  start <- c(mustart, sigmastart, xistart)

  #pass covariate effects on parameters for pp.fit
  n_mu  <- length(attr(model_data$D$mu, "assign"))
  n_sig <- length(attr(model_data$D$sigma, "assign"))
  name  <- c()

  #mul
  if(n_mu == 1){
    mul = NULL
  }else{
    mul = c()
    for (i in 2:n_mu) {
      mu_var <- names(as.data.frame(model_data$D$mu))[i]
      name   <- append(name, mu_var)
      mul    <- append(mul, i-1)
    }
  }

  #sigl
  if(n_sig == 1){
    sigl = NULL
  }else{
    sigl = c()
    for (i in 2:n_sig) {
      sig_var <- names(as.data.frame(model_data$D$sigma))[i]
      name    <- append(name, sig_var)
      sigl     <- append(sigl, i-1)
    }
  }

  y <- model_data$y  #a numeric vector
  if (length(name) == 0) {
    y_thresh <- predict(rq(y ~ 1, p))

  }else{
    index    <- which(colnames(data) == name)
    X        <- data[index]
    y_thresh <- predict(rq(y ~ X, p))

  }

  fit <- ismev::pp.fit(model_data$y, y_thresh, npy, ydat = name,
                       mul = mul, sigl = sigl, shl = NULL,
                       muinit = start[1], siginit = start[2], shinit = start[3],
                       mulink = invmulink, siglink = invsigmalink, shlink = invxilink,
                       show = FALSE)

  # Add information to the returned object
  res <- list()
  res$call         <- Call
  res$data         <- model_data
  res$coefficients <- fit$mle
  #nms <- unlist(lapply(names(model_data$D),
  #                     function(x){
  #                       paste(x, ": ", colnames(model_data$D[[x]]), sep = "")
  #                     }))
  #names(res$coefficients) <- nms
  res$formulae <- mp
  res$invlinks <- list(invmulink = invmulink, invsigmalink = invsigmalink,
                       invxilink = invxilink)
  #res$residuals <-
  res$loglik    <- -fit$nllh
  res$threshold <- y_thresh
  res$cov       <- fit$cov
  res$se        <- fit$se
  res$df        <- length(res$coefficients)
  class(res)    <- c("pp", "evreg")
  res$aic       <- (-2)*res$loglik + 2*(res$df)
  return(res)

}





