#Likelihood-ratio-test-based forward selection

#mu
forward_LRT_mu <- function(fit, alpha = 0.05){

  cov_n <- ncol(eval(fit$call$data)) - 1   #number of covariates in the data

  new_fit <- add1_LRT_mu(fit, alpha)
  re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit

  #Check if the above new_fit is full model
  if(re_n != 1){
    if(length(new_fit$Output_fit$mu) == cov_n){
      return(new_fit)
    }else{
      while (re_n != 1) {
        mu      <- new_fit$Output_fit$mu
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit = add1_LRT_mu(fit, alpha)
      }

      return(new_fit)
    }
  }

}




