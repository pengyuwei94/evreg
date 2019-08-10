#' Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on likelihood-ratio-test,
#' AIC, or p value from Wald test.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param criterion Either based \code{LRT}(Likelihood ratio test),
#'   \code{AIC}, or \code{pvalue}(by Wald test). Default criterion is \code{pvalue}.
#' @param alpha Significance level if criterion equals pvalue or LRT.
#'   Default value is 0.05.
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
#'     \item{added_covariate}{A character vector shows added covariate}
#'     \item{criterion_value}{criterion value for if both input model and output model
#'     are different.}
#' @examples
#' ### Annual Maximum and Minimum Temperature
#'
#' P0 <- gevreg(y = TMX1, data = PORTw[, -1])
#' forward_gevreg(P0)
#'
#' @name forward_gevreg
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname forward_gevreg
#' @export
forward_gevreg <- function(fit, criterion = "pvalue", alpha = 0.05,
                           do_mu = TRUE, do_sigma = FALSE, do_xi = FALSE){
  #1. Check if the input criterion is either LRT, AIC, or pvalue
  if(criterion != "pvalue" && criterion != "AIC" && criterion != "LRT"){
    stop("Criterion should be one of the following: LRT, pvalue and AIC")
  }
  #2. Check if input alpha is valid
  if(alpha > 1 || alpha < 0){
    stop("Significance level alpha should be smaller than 1 and larger than 0")
  }

  do_mu <- eval(do_mu)
  do_sigma <- eval(do_sigma)
  do_xi <- eval(do_xi)
  # Perform forward selection based on specified criterion
  if(criterion == "pvalue"){
    new_fit <- forward_p(fit, alpha, do_mu = do_mu, do_sigma = do_sigma, do_xi = do_xi)
  }else if(criterion == "LRT"){
    new_fit <- forward_LRT(fit, alpha, do_mu = do_mu, do_sigma = do_sigma, do_xi = do_xi)
  }else{
    new_fit <- forward_AIC(fit, do_mu = do_mu, do_sigma = do_sigma, do_xi = do_xi)
  }

  return(new_fit)

}




