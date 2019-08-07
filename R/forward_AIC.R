#' AIC-based Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on AIC.
#'
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
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
#' forward_AIC_mu(P0)
#' forward_AIC_sigma(P0)
#' @name forward_LRT
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname forward_LRT
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
        cov_new   <- length(all.vars(new_fit$Output_fit$mu))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(new_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$mu  <- newer_fit$formulae$mu
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        new_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(new_fit)
}


# ----------------------------- sigma ---------------------------------

#' @rdname forward_LRT
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
        cov_new   <- length(all.vars(new_fit$Output_fit$sigma))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(new_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$sigma <- newer_fit$formulae$sigma
        list$fit   <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        new_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(new_fit)
}


# ----------------------------- xi ---------------------------------

#' @rdname forward_LRT
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
        cov_new   <- length(all.vars(new_fit$Output_fit$xi))
        new_fit   <- newer_fit
      }
      # Make better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call
        newer_fit$AIC <- c(AIC(fit), AIC(newer_fit))
        names(new_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$xi  <- newer_fit$formulae$xi
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list
        # Input_fit
        new_fit$Input_fit <- fit$call

      }


    }
  }else {
    newer_fit <- new_fit
  }

  return(new_fit)
}





