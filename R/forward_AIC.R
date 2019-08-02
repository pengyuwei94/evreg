
#' AIC-based Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on AIC.
#'
#' forward_AIC_mu(fit)
#' forward_AIC_sigma(fit)
#' forward_AIC_xi(fit)
#'
#' @param fit A model of class "gevreg".
#' @details Add details.
#' @return A list which has the following components
#'     \item{Input_fit}{The input object of the class gevreg.}
#'     \item{Note}{A message that will be printed when input fit and output
#'     fit are the same.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class gevreg if output fit is different
#'     from the input fit.}
#'     \item{AIC}{AIC values for both input model and output model if two models
#'     are different.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' forward_AIC_mu(f0)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' library(extRemes)
#' data(PORTw)
#' PORTw$Year <- (PORTw$Year - min(PORTw$Year)) / (max(PORTw$Year) - min(PORTw$Year))
#' P0 <- gevreg(TMX1, data = PORTw)
#'
#' @export
forward_AIC_mu <- function(fit){

  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_AIC_mu(fit)

  # The length of return list from new_fit
  # If this is equal to 3 then we have added a covariate
  # If this is equal to 2 then we haven't
  re_n    <- as.numeric(length(new_fit))

  #Check if the above new_fit is full model
  if(re_n != 2){
    cov_new <- length(all.vars(new_fit$Output_fit$mu))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_LRT_mu().  We stop when either
    # 1. add1_LRT_mu() doesn't add a covariate (re_n = 2), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while re_n != 2 and cov_new < cov_n
    if(cov_new == cov_n){
      return(new_fit)
    }else{
      while (re_n != 2 & cov_new < cov_n) {
        mu      <- new_fit$Output_fit$mu
        model   <- eval(new_fit$Output_fit$fit)
        new_fit <- add1_AIC_mu(model)
        re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit
        cov_new <- length(all.vars(new_fit$Output_fit$mu))
      }
      if(re_n == 2){
        ##Make better output
        # AIC
        new_fit$AIC <- c(AIC(fit), AIC(model))
        names(new_fit$AIC) <- c("Input model", "Output model")
        # Output_fit
        list <- list()
        list$mu  <- mu
        list$fit <- new_fit$Input_fit
        new_fit$Note <- list
        # Input_fit
        new_fit$Input_fit <- fit$call


        # Rename output list
        names(new_fit) <- c("Input_fit", "Output_fit", "AIC")
      }


    }
  }

  return(new_fit)
}
