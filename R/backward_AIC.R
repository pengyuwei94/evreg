#' Backward Elimination on GEV Parameter based on AIC
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with backward direction based on AIC.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param do_mu do backward selection on mu if \code{do_mu} equals TRUE. Default is TRUE.
#' @param do_sigma do backward selection on sigma if \code{do_sigma} equals TRUE. Default is FALSE.
#' @param do_xi do backward selection on xi if \code{do_xi} equals TRUE. Default is FALSE.
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
#'     \item{AIC}{AIC values for both input model and output model if two models
#'     are different.}
#' @examples
#'
#' ### Oxford and Worthing annual maximum temperatures
#'
#' ow$year <- (ow$year - 1901) / (1980 - 1901)
#' ow1 <- gevreg(y = temp, data = ow[-3], mu = ~loc + year, sigma = ~loc,
#' xi = ~loc, sigmalink = identity)
#' backward_AIC_mu(ow1)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P6 <- gevreg(TMX1, data = PORTw[,-1], mu = ~MTMAX + AOindex + STDTMAX + STDMIN + MDTR)
#' P7 <- gevreg(TMX1, data = PORTw[,-1], sigma = ~MTMAX + STDTMAX + STDMIN + MDTR)
#' P8 <- gevreg(TMX1, data = PORTw[,-1], xi = ~MTMAX + STDTMAX + STDMIN + MDTR)
#' backward_AIC_mu(P6)
#' backward_AIC_sigma(P7)
#' backward_AIC_xi(P8)
#'
#' @name backward_AIC
NULL
## NULL


#' @rdname backward_AIC
#' @export
backward_AIC <- function(fit, do_mu = TRUE, do_sigma = FALSE, do_xi = FALSE){
  #1. If only performing backward selection on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == FALSE){
    AIC_mu  <- backward_AIC_mu(fit)
    new_fit <- AIC_mu
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #2. If performing backward selection first on sigma, then mu
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == FALSE){
    AIC_sigma <- backward_AIC_sigma(fit)
    AIC_mu    <- backward_AIC_mu(AIC_sigma)
    new_fit   <- AIC_mu
    new_fit$dropped_covariate <- append(AIC_sigma$dropped_covariate, AIC_mu$dropped_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #3. If performing backward selection first on xi, second on sigma, then on mu
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == TRUE){
    AIC_xi    <- backward_AIC_xi(fit)
    AIC_sigma <- backward_AIC_sigma(AIC_mu)
    AIC_mu    <- backward_AIC_mu(AIC_sigma)
    new_fit   <- AIC_mu
    new_fit$dropped_covariate <- append(AIC_xi$dropped_covariate,
                                      AIC_sigma$dropped_covariate)
    new_fit$dropped_covariate <- append(new_fit$dropped_covariate,
                                      AIC_mu$dropped_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #4. If performing backward selection first on xi, then on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == TRUE){
    AIC_xi  <- backward_AIC_xi(fit)
    AIC_mu  <- backward_AIC_mu(AIC_xi)
    new_fit <- AIC_mu
    new_fit$dropped_covariate <- append(AIC_xi$dropped_covariate, AIC_mu$dropped_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #5. If performing backward selection first on xi, then on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == TRUE){
    AIC_xi  <- backward_AIC_xi(fit)
    AIC_sigma  <- backward_AIC_mu(AIC_xi)
    new_fit <- AIC_sigma
    new_fit$dropped_covariate <- append(AIC_xi$dropped_covariate, AIC_sigma$dropped_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #6. If performing backward selection only on xi
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == TRUE){
    AIC_xi  <- backward_AIC_xi(fit)
    new_fit <- AIC_xi
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #7. If performing backward selection only on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == FALSE){
    AIC_sigma  <- backward_AIC_sigma(fit)
    new_fit <- AIC_sigma
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #8. If performing no backward selection on any parameters
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == FALSE){
    new_fit <- fit
  }
  return(new_fit)
}


# ----------------------------- mu ---------------------------------

#' @rdname backward_AIC
#' @export
backward_AIC_mu <- function(fit){

  new_fit <- drop1_AIC_mu(fit)
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$mu)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to dropp more variables, one at a time, using
    # drop1_AIC_mu().  We stop when either
    # 1. drop1_AIC_mu() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_AIC_mu(new_fit)
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

        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")

      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}


# ----------------------------- sigma ---------------------------------

#' @rdname backward_AIC
#' @export
backward_AIC_sigma <- function(fit){

  new_fit <- drop1_AIC_sigma(fit)
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$sigma)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to dropp more variables, one at a time, using
    # drop1_AIC_sigma().  We stop when either
    # 1. drop1_AIC_sigma() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_AIC_sigma(new_fit)
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

        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")

      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}



# ----------------------------- xi ---------------------------------

#' @rdname backward_AIC
#' @export
backward_AIC_xi <- function(fit){

  new_fit <- drop1_AIC_xi(fit)
  # Store dropped_covariate for later use
  dropped <- new_fit$dropped_covariate
  #Check if the above new_fit is null model
  if(new_fit$Note == "covariate dropped"){
    cov_new <- ncol(new_fit$data$D$xi)
    # If we have dropped all the covariates then we stop
    # Otherwise, we try to dropp more variables, one at a time, using
    # drop1_AIC_xi().  We stop when either
    # 1. drop1_AIC_xi() doesn't drop a covariate (new_fit$Note != "covariate dropped"), or
    # 2. we have dropped all the covariates (cov_new == 1)
    # Therefore, we continue to loop while new_fit$Note == "covariate dropped" and cov_new != 1
    if(cov_new == 1){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (new_fit$Note == "covariate dropped" & cov_new != 1) {
        newer_fit <- drop1_AIC_xi(new_fit)
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

        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")

      }
    }
  }else{
    newer_fit <- new_fit
  }

  return(newer_fit)
}

