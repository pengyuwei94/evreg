#' Backward Elimination on GEV Parameter based on Likelihood-ratio-test
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with backward direction based on Likelihood-ratio-test.
#'
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if covariates have been dropped or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{dropped_covariate}{A character vector shows dropped covariates}
#'     \item{pvalue}{A data frame that contains p values with five decimal
#'     places of the Likelihood-ratio-test.}
#' @examples
#' ### Oxford and Worthing annual maximum temperatures
#'
#' ow$year <- (pjn$year - 1901) / (1980 - 1901)
#' ow1 <- gevreg(y = temp, data = ow[-3], mu = ~loc + year, sigma = ~loc,
#' xi = ~loc, sigmalink = identity)
#' backward_LRT_xi(ow1)
#'
#'
#' #' ### Annual Maximum and Minimum Temperature
#'
#' P3 <- gevreg(y = TMX1, data = PORTw[, -1], mu = ~MTMAX + STDTMAX + STDMIN)
#' P5 <- gevreg(y = TMX1, data = PORTw[, -1], xi = ~MTMAX + AOindex)
#' backward_LRT_mu(P3)
#' backward_LRT_xi(P5)
#'
#' @name backward_LRT
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname backward_LRT
#' @export
backward_LRT_mu <- function(fit, alpha = 0.05){

  new_fit <- drop1_LRT_mu(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$mu)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to drop more variables, one at a time, using
    # drop1_LRT_mu().  We stop when either
    # 1. drop1_LRT_mu() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_LRT_mu(new_fit, alpha)
        p_table <- rbind(p_table, newer_fit$pvalue)
        dropped <- append(dropped, newer_fit$dropped_covariate)
        cov_new <- ncol(newer_fit$data$D$mu)
        new_fit <- newer_fit
      }
      # Make better output
      if(new_fit$Note != "covariate dropped"){
        newer_fit$dropped_covariate <- dropped
        newer_fit$Note <- "covariate dropped"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$mu  <- newer_fit$formulae$mu
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)
      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}


# ----------------------------- sigma ---------------------------------

#' @rdname backward_LRT
#' @export
backward_LRT_sigma <- function(fit, alpha = 0.05){

  new_fit <- drop1_LRT_sigma(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$sigma)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to drop more variables, one at a time, using
    # drop1_LRT_sigma().  We stop when either
    # 1. drop1_LRT_sigma() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_LRT_sigma(new_fit, alpha)
        p_table <- rbind(p_table, newer_fit$pvalue)
        dropped <- append(dropped, newer_fit$dropped_covariate)
        cov_new <- ncol(newer_fit$data$D$sigma)
        new_fit <- newer_fit
      }
      # Make better output
      if(new_fit$Note != "covariate dropped"){
        newer_fit$dropped_covariate <- dropped
        newer_fit$Note <- "covariate dropped"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$sigma <- newer_fit$formulae$sigma
        list$fit   <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)
      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}


# ----------------------------- xi ---------------------------------

#' @rdname backward_LRT
#' @export
backward_LRT_xi <- function(fit, alpha = 0.05){

  new_fit <- drop1_LRT_xi(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$xi)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to drop more variables, one at a time, using
    # drop1_LRT_xi().  We stop when either
    # 1. drop1_LRT_xi() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_LRT_xi(new_fit, alpha)
        p_table <- rbind(p_table, newer_fit$pvalue)
        dropped <- append(dropped, newer_fit$dropped_covariate)
        cov_new <- ncol(newer_fit$data$D$xi)
        new_fit <- newer_fit
      }
      # Make better output
      if(new_fit$Note != "covariate dropped"){
        newer_fit$dropped_covariate <- dropped
        newer_fit$Note <- "covariate dropped"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$xi  <- newer_fit$formulae$xi
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)
      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}


