#' Likelihood-ratio-test-based Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on likelihood-ratio-test.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05.
#' @param do_mu do forward selection on mu if \code{do_mu} equals TRUE. Default is TRUE.
#' @param do_sigma do forward selection on sigma if \code{do_sigma} equals TRUE. Default is FALSE.
#' @param do_xi do forward selection on xi if \code{do_xi} equals TRUE. Default is FALSE.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if covariates have been added or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{added_covariate}{A character vector shows added covariate}
#'     \item{pvalue}{A data frame that contains p values with five decimal
#'     places of the Likelihood-ratio-test.}
#' @examples
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle[, -1])
#' forward_LRT_mu(f0)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P0 <- gevreg(y = TMX1, data = PORTw[, -1])
#' forward_LRT(P0)

#'
#' @name forward_LRT
NULL
## NULL


#' @rdname forward_LRT
#' @export
forward_LRT <- function(fit, alpha = 0.05,
                        do_mu = TRUE, do_sigma = FALSE, do_xi = FALSE){
  #1. If only performing forward selection on mu
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == FALSE){
    LRT_mu  <- forward_LRT_mu(fit)
    new_fit <- LRT_mu
  }
  #2. If performing forward selection first on mu, then sigma
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == FALSE){
    LRT_mu    <- forward_LRT_mu(fit)
    LRT_sigma <- forward_LRT_sigma(LRT_mu)
    new_fit   <- LRT_sigma
   # Make better output for criterion value
    new_fit$added_covariate <- append(LRT_mu$added_covariate, LRT_sigma$added_covariate)
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note == "covariate added"){
     new_fit$pvalue              <- rbind(LRT_mu$pvalue, LRT_sigma$pvalue)
    }
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note != "covariate added"){
     new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note != "covariate added" && LRT_sigma$Note == "covariate added"){
     new_fit$pvalue              <- LRT_sigma$pvalue
    }
  }
  #3. If performing forward selection first on mu, second on sigma, then on xi
  if(do_mu == TRUE && do_sigma == TRUE && do_xi == TRUE){
    LRT_mu    <- forward_LRT_mu(fit)
    LRT_sigma <- forward_LRT_sigma(LRT_mu)
    LRT_xi    <- forward_LRT_xi(LRT_sigma)
    new_fit   <- LRT_xi
    # Make better output for criterion value
    new_fit$added_covariate <- append(LRT_mu$added_covariate,
                                      LRT_sigma$added_covariate)
    new_fit$added_covariate <- append(new_fit$added_covariate,
                                      LRT_xi$added_covariate)
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note == "covariate added" && LRT_xi$Note == "covariate added"){
     new_fit$pvalue              <- rbind(LRT_mu$pvalue, LRT_sigma$pvalue,  LRT_xi$pvalue)
    }
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note != "covariate added" && LRT_xi$Note != "covariate added"){
     new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note == "covariate added" && LRT_xi$Note != "covariate added"){
     new_fit$pvalue              <- rbind(LRT_mu$pvalue, LRT_sigma$pvalue)
    }
    if(LRT_mu$Note == "covariate added" && LRT_sigma$Note != "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- rbind(LRT_mu$pvalue, LRT_xi$pvalue)
    }
    if(LRT_mu$Note != "covariate added" && LRT_sigma$Note != "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
    if(LRT_mu$Note != "covariate added" && LRT_sigma$Note == "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- rbind(LRT_sigma$pvalue, LRT_xi$pvalue)
    }
    if(LRT_mu$Note != "covariate added" && LRT_sigma$Note == "covariate added" && LRT_xi$Note != "covariate added"){
      new_fit$pvalue              <- LRT_sigma$pvalue
    }
  }
  #4. If performing forward selection first on mu, then on xi
  if(do_mu == TRUE && do_sigma == FALSE && do_xi == TRUE){
    LRT_mu  <- forward_LRT_mu(fit)
    LRT_xi  <- forward_LRT_xi(LRT_mu)
    new_fit <- LRT_xi
    # Make better output for criterion value
    new_fit$added_covariate <- append(LRT_mu$added_covariate, LRT_xi$added_covariate)
    if(LRT_mu$Note == "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- rbind(LRT_mu$pvalue, LRT_xi$pvalue)
    }
    if(LRT_mu$Note == "covariate added" && LRT_xi$Note != "covariate added"){
      new_fit$pvalue              <- LRT_mu$pvalue
    }
    if(LRT_mu$Note != "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
  }
  #5. If performing forward selection first on sigma, then on xi
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == TRUE){
    LRT_sigma  <- forward_LRT_mu(fit)
    LRT_xi  <- forward_LRT_xi(LRT_sigma)
    new_fit <- LRT_xi
    # Make better output for criterion value
    new_fit$added_covariate <- append(LRT_sigma$added_covariate, LRT_xi$added_covariate)
    if(LRT_sigma$Note == "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- rbind(LRT_sigma$pvalue, LRT_xi$pvalue)
    }
    if(LRT_sigma$Note == "covariate added" && LRT_xi$Note != "covariate added"){
      new_fit$pvalue              <- LRT_sigma$pvalue
    }
    if(LRT_sigma$Note != "covariate added" && LRT_xi$Note == "covariate added"){
      new_fit$pvalue              <- LRT_xi$pvalue
    }
  }
  #6. If performing forward selection only on xi
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == TRUE){
    LRT_xi  <- forward_LRT_xi(fit)
    new_fit <- LRT_xi
  }
  #7. If performing forward selection only on sigma
  if(do_mu == FALSE && do_sigma == TRUE && do_xi == FALSE){
    LRT_sigma  <- forward_LRT_sigma(fit)
    new_fit <- LRT_sigma
  }
  #8. If performing no forward selection on any parameters
  if(do_mu == FALSE && do_sigma == FALSE && do_xi == FALSE){
    new_fit <- fit
  }
  return(new_fit)
}


# ----------------------------- mu ---------------------------------

#' @rdname forward_LRT
#' @export
forward_LRT_mu <- function(fit, alpha = 0.05){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_LRT_mu(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store added_covariate for later use
  added   <- new_fit$added_covariate
  # Check if the above new_fit is full model
  if (new_fit$Note == "covariate added") {
    cov_new <- length(all.vars(new_fit$formulae$mu))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_LRT_mu().  We stop when either
    # 1. add1_LRT_mu() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if (cov_new == cov_n) {
      newer_fit <- new_fit
    } else {
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_LRT_mu(new_fit, alpha)
        p_table   <- rbind(p_table, newer_fit$pvalue)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$formulae$mu))
        new_fit   <- newer_fit
      }
      # Make a better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$mu  <- newer_fit$formulae$mu
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)

      }

    }
  } else {
    newer_fit <- new_fit
  }
  return(newer_fit)
}



# ----------------------------- sigma ---------------------------------

#' @rdname forward_LRT
#' @export
forward_LRT_sigma <- function(fit, alpha = 0.05){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_LRT_sigma(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store added_covariate for later use
  added   <- new_fit$added_covariate
  # Check if the above new_fit is full model
  if (new_fit$Note == "covariate added") {
    cov_new <- length(all.vars(new_fit$formulae$sigma))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_LRT_sigma().  We stop when either
    # 1. add1_LRT_sigma() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if (cov_new == cov_n) {
      newer_fit <- new_fit
    } else {
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_LRT_sigma(new_fit, alpha)
        p_table   <- rbind(p_table, newer_fit$pvalue)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$formulae$sigma))
        new_fit   <- newer_fit
      }
      # Make a better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$sigma <- newer_fit$formulae$sigma
        list$fit   <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)

      }

    }
  } else {
    newer_fit <- new_fit
  }
  return(newer_fit)
}


# ----------------------------- xi ---------------------------------

#' @rdname forward_LRT
#' @export
forward_LRT_xi <- function(fit, alpha = 0.05){
  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_LRT_xi(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue
  # Store added_covariate for later use
  added   <- new_fit$added_covariate
  # Check if the above new_fit is full model
  if (new_fit$Note == "covariate added") {
    cov_new <- length(all.vars(new_fit$formulae$sigma))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_LRT_xi().  We stop when either
    # 1. add1_LRT_xi() doesn't add a covariate (new_fit$Note != "covariate added"), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while new_fit$Note == "covariate added" and cov_new < cov_n
    if (cov_new == cov_n) {
      newer_fit <- new_fit
    } else {
      newer_fit <- new_fit
      while (newer_fit$Note == "covariate added" & cov_new < cov_n) {
        newer_fit <- add1_LRT_xi(new_fit, alpha)
        p_table   <- rbind(p_table, newer_fit$pvalue)
        added     <- append(added, newer_fit$added_covariate)
        cov_new   <- length(all.vars(newer_fit$formulae$xi))
        new_fit   <- newer_fit
      }
      # Make a better output
      if(newer_fit$Note != "covariate added"){
        newer_fit$added_covariate <- added
        newer_fit$Note <- "covariate added"
        newer_fit$Input_fit <- fit$call

        list <- list()
        list$xi  <- newer_fit$formulae$xi
        list$fit <- newer_fit$call
        newer_fit$Output_fit <- list

        newer_fit$pvalue <- head(p_table,-1)

      }

    }
  } else {
    newer_fit <- new_fit
  }
  return(newer_fit)
}



