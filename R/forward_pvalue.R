#' Forward Selection on GEV Parameter using individual p value
#'
#' Significance controlled variable selection selects variables in either
#' mu, sigma, and xi with forward direction based on individual p value.
#'
#' forward_p_mu(fit, alpha = 0.05)
#' forward_p_sigma(fit, alpha = 0.05)
#' forward_p_xi(fit, alpha = 0.05)
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
#'     places of the added covariates in order if there are.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' forward_p_mu(f0)
#'
#' @export
forward_p_mu <- function(fit, alpha = 0.05){

  # Number of covariates in the data
  cov_n   <- ncol(eval(fit$call$data)) - 1
  new_fit <- add1_p_mu(fit, alpha)

  # Store p table for later use
  p_table <- new_fit$pvalue


  # The length of return list from new_fit
  # If this is equal to 3 then we have added a covariate
  # If this is equal to 2 then we haven't
  re_n    <- as.numeric(length(new_fit))

  #Check if the above new_fit is full model
  if(re_n != 2){
    cov_new <- length(all.vars(new_fit$Output_fit$mu))
    # If we have added all the covariates then we stop
    # Otherwise, we try adding more variables, one at a time, using
    # add1_p_mu().  We stop when either
    # 1. add1_p_mu() doesn't add a covariate (re_n = 2), or
    # 2. we have added all the covariates (cov_new = cov_n)
    # Therefore, we continue to loop while re_n != 2 and cov_new < cov_n
    if(cov_new == cov_n){
      return(new_fit)
    }else{
      while (re_n != 2 & cov_new < cov_n) {
        mu      <- new_fit$Output_fit$mu
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit <- add1_p_mu(fit, alpha)
        # Improve the p value data frame from add1_p_mu()
        new_p   <- new_fit$pvalue
        p_table <- rbind(p_table, new_p)

        re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit
        cov_new <- length(all.vars(new_fit$Output_fit$mu))
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



