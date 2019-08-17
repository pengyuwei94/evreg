#' AIC-based Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on AIC.
#'
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param do_mu do forward selection on mu if \code{do_mu} equals TRUE. Default is TRUE.
#' @param do_sigma do forward selection on sigma if \code{do_sigma} equals TRUE. Default is FALSE.
#' @param do_xi do forward selection on xi if \code{do_xi} equals TRUE. Default is FALSE.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if a covariate has been added or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{added_covariate}{A character vector shows added covariates}
#'     \item{AIC}{AIC values for both input model and output model if two models
#'     are different.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle[,-1])
#' forward_AIC_mu(f0)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P0 <- gevreg(y = TMX1, data = PORTw[, -1])
#' forward_AIC(P0)
#' @name forward_AIC
NULL
## NULL


#' @rdname forward_AIC
#' @export
forward_AIC <- function(fit, do_mu = TRUE, do_sigma = FALSE, do_xi = FALSE){
  #1. If only performing forward selection on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == FALSE){
    AIC_mu  <- forward_AIC_mu(fit)
    new_fit <- AIC_mu
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #2. If performing forward selection first on mu, then sigma
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == FALSE){
    AIC_mu    <- forward_AIC_mu(fit)
    AIC_sigma <- forward_AIC_sigma(AIC_mu)
    new_fit   <- AIC_sigma
    new_fit$added_covariate <- append(AIC_mu$added_covariate, AIC_sigma$added_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #3. If performing forward selection first on mu, second on sigma, then on xi
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == TRUE){
    AIC_mu    <- forward_AIC_mu(fit)
    AIC_sigma <- forward_AIC_sigma(AIC_mu)
    AIC_xi    <- forward_AIC_xi(AIC_sigma)
    new_fit   <- AIC_xi
    new_fit$added_covariate <- append(AIC_mu$added_covariate,
                                      AIC_sigma$added_covariate,
                                      AIC_xi$added_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #4. If performing forward selection first on mu, then on xi
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == TRUE){
    AIC_mu  <- forward_AIC_mu(fit)
    AIC_xi  <- forward_AIC_xi(AIC_mu)
    new_fit <- AIC_xi
    new_fit$added_covariate <- append(AIC_mu$added_covariate, AIC_xi$added_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #5. If performing forward selection first on sigma, then on xi
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == TRUE){
    AIC_sigma  <- forward_AIC_mu(fit)
    AIC_xi  <- forward_AIC_xi(AIC_sigma)
    new_fit <- AIC_xi
    new_fit$added_covariate <- append(AIC_sigma$added_covariate, AIC_xi$added_covariate)
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #6. If performing forward selection only on xi
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == TRUE){
    AIC_xi  <- forward_AIC_xi(fit)
    new_fit <- AIC_xi
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #7. If performing forward selection only on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == FALSE){
    AIC_sigma  <- forward_AIC_sigma(fit)
    new_fit <- AIC_sigma
    new_fit$AIC <- c(AIC(fit), AIC(new_fit))
    names(new_fit$AIC) <- c("Input model", "Output model")
  }
  #8. If performing no forward selection on any parameters
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == FALSE){
    new_fit <- fit
  }
  return(new_fit)
}


# ----------------------------- mu ---------------------------------

#' @rdname forward_AIC
#' @export
forward_AIC_mu <- function(fit){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_AIC_mu(fit)
  # Store added_covariate for later use
  added   <- new_fit$added_covariate

  #Check if the above new_fit is full model
  if(new_fit$Note == "covariate added"){
    cov_new <- length(all.vars(new_fit$formulae$mu))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_AIC_mu().  We stop when either
    # 1. add1_AIC_mu() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if(cov_new == cov_n){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_AIC_mu(new_fit)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$Output_fit$mu))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$mu  <- newer_fit$formulae$mu
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        newer_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(newer_fit)
}


# ----------------------------- sigma ---------------------------------

#' @rdname forward_AIC
#' @export
forward_AIC_sigma <- function(fit){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_AIC_sigma(fit)
  # Store added_covariate for later use
  added   <- new_fit$added_covariate

  #Check if the above new_fit is full model
  if(new_fit$Note == "covariate added"){
    cov_new <- length(all.vars(new_fit$formulae$sigma))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_AIC_sigma().  We stop when either
    # 1. add1_AIC_sigma() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if(cov_new == cov_n){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_AIC_sigma(new_fit)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$Output_fit$sigma))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$sigma <- newer_fit$formulae$sigma
        list$fit   <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        newer_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(newer_fit)
}


# ----------------------------- xi ---------------------------------

#' @rdname forward_AIC
#' @export
forward_AIC_xi <- function(fit){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_AIC_xi(fit)
  # Store added_covariate for later use
  added   <- new_fit$added_covariate

  #Check if the above new_fit is full model
  if(new_fit$Note == "covariate added"){
    cov_new <- length(all.vars(new_fit$formulae$xi))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_AIC_xi().  We stop when either
    # 1. add1_AIC_xi() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if(cov_new == cov_n){
      newer_fit <- new_fit
    }else{
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_AIC_xi(new_fit)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$Output_fit$xi))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(newer_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$xi  <- newer_fit$formulae$xi
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        newer_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(newer_fit)
}





