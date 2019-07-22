#AIC-based forward selection

#mu
forward_AIC_mu <- function(fit){
  cov_n <- ncol(eval(fit$call$data)) - 1   #number of covariates in the data

  new_fit <- add1_AIC_mu(fit)
  re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit

  #Check if the above new_fit is full model
  if(re_n != 1){
    if(length(new_fit$Output_fit$mu) == cov_n){
      return(new_fit)
    }else{
      while (re_n != 1) {
        mu      <- new_fit$Output_fit$mu
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit <- add1_AIC_mu(fit)
      }

      return(new_fit)
    }
  }


}
