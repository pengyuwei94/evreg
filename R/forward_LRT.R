
#' Likelihood-ratio-test-based Forward Selection on GEV Parameter
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on likelihood-ratio-test.
#'
#' forward_LRT_mu(fit, alpha = 0.05)
#' forward_LRT_sigma(fit, alpha = 0.05)
#' forward_LRT_xi(fit, alpha = 0.05)
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
#'     \item{pvalue}{A p value based on a likelihood-ratio-test if the input fit
#'     and output fit are different.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' forward_LRT_mu(f0)
#'
#' @export
forward_LRT_mu <- function(fit, alpha = 0.05){

  # Number of covariates in the data
  cov_n <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_LRT_mu(fit, alpha)
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
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit <- add1_LRT_mu(fit, alpha)
        re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit
        cov_new <- length(all.vars(new_fit$Output_fit$mu))
      }

    }
  }

  return(new_fit)
}




