#' Backward Elimination on GEV Parameter based on Likelihood-ratio-test
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with backward direction based on Likelihood-ratio-test.
#'
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05.
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

#' @name backward_LRT
NULL
## NULL


#' @rdname backward_LRT
#' @export
backward_LRT <- function(fit, alpha = 0.05,
                        do_mu = TRUE, do_sigma = FALSE, do_xi = FALSE){
  #1. If only performing backward selection on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == FALSE){
    LRT_mu  <- backward_LRT_mu(fit)
    new_fit <- LRT_mu
  }
  #2. If performing backward selection first on sigma, then mu
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == FALSE){
    LRT_sigma <- backward_LRT_sigma(fit)
    LRT_mu    <- backward_LRT_mu(LRT_sigma)
    new_fit   <- LRT_mu
    # Make better output for criterion value
    new_fit$dropped_covariate <- append(LRT_sigma$dropped_covariate, LRT_mu$dropped_covariate)
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_sigma$pvalue, LRT_mu$pvalue)
    }
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note != "covariate dropped"){
      new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note != "covariate dropped" && LRT_sigma$Note == "covariate dropped"){
      new_fit$pvalue              <- LRT_sigma$pvalue
    }
  }
  #3. If performing backward selection first on xi, second on sigma, then on mu
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == TRUE){
    LRT_xi    <- backward_LRT_xi(fit)
    LRT_sigma <- backward_LRT_sigma(LRT_xi)
    LRT_mu    <- backward_LRT_mu(LRT_sigma)
    new_fit   <- LRT_mu
    # Make better output for criterion value
    new_fit$dropped_covariate <- append(LRT_xi$dropped_covariate,
                                      LRT_sigma$dropped_covariate,
                                      LRT_mu$dropped_covariate)
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note == "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_xi$pvalue, LRT_sigma$pvalue,  LRT_mu$pvalue)
    }
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note != "covariate dropped" && LRT_xi$Note != "covariate dropped"){
      new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note == "covariate dropped" && LRT_xi$Note != "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_sigma$pvalue, LRT_mu$pvalue)
    }
    if(LRT_mu$Note == "covariate dropped" && LRT_sigma$Note != "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_xi$pvalue, LRT_mu$pvalue)
    }
    if(LRT_mu$Note != "covariate dropped" && LRT_sigma$Note != "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
    if(LRT_mu$Note != "covariate dropped" && LRT_sigma$Note == "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_xi$pvalue, LRT_sigma$pvalue)
    }
    if(LRT_mu$Note != "covariate dropped" && LRT_sigma$Note == "covariate dropped" && LRT_xi$Note != "covariate dropped"){
      new_fit$pvalue              <- LRT_sigma$pvalue
    }
  }
  #4. If performing backward selection first on xi, then on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == TRUE){
    LRT_xi  <- backward_LRT_xi(fit)
    LRT_mu  <- backward_LRT_mu(LRT_xi)
    new_fit <- LRT_mu
    # Make better output for criterion value
    new_fit$dropped_covariate <- append(LRT_xi$dropped_covariate, LRT_mu$dropped_covariate)
    if(LRT_mu$Note == "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_xi$pvalue, LRT_mu$pvalue)
    }
    if(LRT_mu$Note == "covariate dropped" && LRT_xi$Note != "covariate dropped"){
      new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note != "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
  }
  #5. If performing backward selection first on xi, then on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == TRUE){
    LRT_xi  <- backward_LRT_xi(fit)
    LRT_sigma  <- backward_LRT_mu(LRT_xi)
    new_fit <- LRT_sigma
    # Make better output for criterion value
    new_fit$dropped_covariate <- append(LRT_xi$dropped_covariate, LRT_sigma$dropped_covariate)
    if(LRT_sigma$Note == "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- rbind(LRT_xi$pvalue, LRT_sigma$pvalue)
    }
    if(LRT_sigma$Note == "covariate dropped" && LRT_xi$Note != "covariate dropped"){
      new_fit$pvalue              <- LRT_sigma$pvalue
    }
    if(LRT_sigma$Note != "covariate dropped" && LRT_xi$Note == "covariate dropped"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
  }
  #6. If performing backward selection only on xi
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == TRUE){
    LRT_xi  <- backward_LRT_xi(fit)
    new_fit <- LRT_xi
  }
  #7. If performing backward selection only on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == FALSE){
    LRT_sigma  <- backward_LRT_sigma(fit)
    new_fit <- LRT_sigma
  }
  #8. If performing no backward selection on any parameters
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == FALSE){
    new_fit <- fit
  }
  return(new_fit)
}


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


