#' Backward Elimination on GEV Parameter using individual p value
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with backward direction based on individual p value.
#'
#' backward_p_mu(fit, alpha = 0.05)
#' backward_p_sigma(fit, alpha = 0.05)
#' backward_p_xi(fit, alpha = 0.05)
#'
#' @param fit A model of class "gevreg".
#' @param alpha Significance level. Default value is 0.05..
#' @details Add details.
#' @return A list which has the following components
#'     \item{Input_fit}{The input object of the class gevreg.}
#'     \item{Note}{A message that will be printed when input fit and output
#'     fit are the same.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class gevreg if output fit is different
#'     from the input fit.}
#'     \item{pvalue}{A data frame that contains p values with five decimal
#'     places of the dropped covariates in order if there are.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f3 <- gevreg(y = SeaLevel, data = evreg::fremantle[-1], mu = ~Year01 + SOI)
#' backward_p_mu(f3)
#'
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
#' backward_p_mu(P6)
#'
#' @export
backward_p_mu <- function(fit, alpha = 0.05){

  new_fit <- drop1_p_mu(fit, alpha)
  # Store p table for later use
  p_table <- new_fit$pvalue

  # The length of return list from new_fit
  # If this is equal to 3 then we have added a covariate
  # If this is equal to 2 then we haven't
  re_n    <- as.numeric(length(new_fit))

  #Check if the above new_fit is null model
  if(re_n != 2){
    var <- as.character(new_fit$Output_fit$mu)[2]

    # If we have dropped all the covariates then we stop
    # Otherwise, we try to dropp more variables, one at a time, using
    # drop1_p_mu().  We stop when either
    # 1. drop1_p_mu() doesn't drop a covariate (re_n = 2), or
    # 2. we have dropped all the covariates (var == "1")
    # Therefore, we continue to loop while re_n != 2 and var != "1"
    if(var == "1"){
      return(new_fit)
    }else{
      while (re_n != 2 & var != "1") {
        mu      <- new_fit$Output_fit$mu
        model   <- eval(new_fit$Output_fit$fit)
        new_fit <- drop1_p_mu(model, alpha)
        # Improve the p value data frame from drop1_p_mu()
        new_p   <- new_fit$pvalue
        p_table <- rbind(p_table, new_p)

        re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit
        var     <- as.character(new_fit$Output_fit$mu)[2]
      }
      if(re_n == 2){
        ##Make better output
        # Output_fit
        list <- list()
        list$mu  <- mu
        list$fit <- new_fit$Input_fit
        new_fit$Note <- list
        # Input_fit
        new_fit$Input_fit <- fit$call
        # Update p table
        new_fit$pvalue <- p_table

        # Rename output list
        names(new_fit) <- c("Input_fit", "Output_fit", "pvalue")
      }
    }


  }

  return(new_fit)
}


