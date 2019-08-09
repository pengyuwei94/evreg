#' Add one possible covariate on GEV parameter
#'
#' Add a single term to either mu, sigma, and xi based on criterion
#' Likelihood ratio test, AIC, or p value from Wald test
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param parameter A specified parameter that needs to add one covariate.
#'   Equals \code{"mu"} by default.
#' @param criterion Either based \code{LRT}(Likelihood ratio test),
#'   \code{AIC}, or \code{pvalue}(by Wald test). Default criterion is \code{pvalue}.
#' @param alpha Significance level if criterion equals pvalue or LRT.
#'   Default value is 0.05.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if a covariate has been added or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{added_covariate}{A character vector shows added covariate}
#'     \item{criterion_value}{criterion value for if both input model and output model
#'     are different.}
#' @examples
#' ### Annual Maximum and Minimum Temperature
#'
#' P0 <- gevreg(y = TMX1, data = PORTw[, -1])
#' add1_gevreg(P0)
#'
#' @name add1_gevreg
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname add1_gevreg
#' @export
add1_gevreg <- function(fit, parameter = "mu", alpha = 0.05, criterion = "pvalue"){
  #1. Check if the input criterion is either LRT, AIC, or pvalue
  if(criterion != "pvalue" && criterion != "AIC" && criterion != "LRT"){
    stop("Criterion should be one of the following: LRT, pvalue and AIC")
  }
  #2. Check if input alpha is valid
  if(alpha > 1 || alpha < 0){
    stop("Significance level alpha should be smaller than 1 and larger than 0")
  }
  #3. Check if the input parameter is either mu, sigma, or xi
  if(parameter != "mu" && parameter != "sigma" && parameter != "xi"){
    stop("Parameter should be one of the following: mu, sigma, or xi")
  }
  # Specify criterion
  if(criterion == "pvalue"){
    #Specify parameter
    if(parameter == "mu"){
      new_fit <- add1_p_mu(fit)
    }else if(parameter == "sigma"){
      new_fit <- add1_p_sigma(fit)
    }else{
      new_fit <- add1_p_xi(fit)
    }
  }else if(criterion == "AIC"){
    #Specify parameter
    if(parameter == "mu"){
      new_fit <- add1_AIC_mu(fit)
    }else if(parameter == "sigma"){
      new_fit <- add1_AIC_sigma(fit)
    }else{
      new_fit <- add1_AIC_xi(fit)
    }
  }else{
    #Specify parameter
    if(parameter == "mu"){
      new_fit <- add1_LRT_mu(fit)
    }else if(parameter == "sigma"){
      new_fit <- add1_LRT_sigma(fit)
    }else{
      new_fit <- add1_LRT_xi(fit)
    }
  }

  return(new_fit)
}



