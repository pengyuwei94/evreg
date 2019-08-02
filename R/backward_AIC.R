#' Backward Elimination on GEV Parameter based on AIC
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with backward direction based on AIC.
#'
#' backward_AIC_mu(fit)
#' backward_AIC_sigma(fit)
#' backward_AIC_xi(fit)
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
#' library(extRemes)
#' data(PORTw)
#' PORTw$Year <- (PORTw$Year - min(PORTw$Year)) / (max(PORTw$Year) - min(PORTw$Year))
#' P6 <- gevreg(TMX1, data = PORTw, mu = ~MTMAX + AOindex + Year + STDTMAX + STDMIN + MDTR)
#' backward_AIC_mu(P6)
#'
#' @export
backward_AIC_mu <- function(fit){

  new_fit <- drop1_AIC_mu(fit)
  # The length of return list from new_fit
  # If this is equal to 3 then we have added a covariate
  # If this is equal to 2 then we haven't
  re_n    <- as.numeric(length(new_fit))

  #Check if the above new_fit is null model
  if(re_n != 2){
    var <- as.character(new_fit$Output_fit$mu)[2]

    # If we have dropped all the covariates then we stop
    # Otherwise, we try to dropp more variables, one at a time, using
    # drop1_AIC_mu().  We stop when either
    # 1. drop1_AIC_mu() doesn't drop a covariate (re_n = 2), or
    # 2. we have dropped all the covariates (var == "1")
    # Therefore, we continue to loop while re_n != 2 and var != "1"
    if(var == "1"){
      return(new_fit)
    }else{
      while (re_n != 2 & var != "1") {
        mu      <- new_fit$Output_fit$mu
        model   <- eval(new_fit$Output_fit$fit)
        new_fit <- drop1_AIC_mu(model)

        re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit
        var     <- as.character(new_fit$Output_fit$mu)[2]
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
