#Likelihood-ratio-test-based backward elimination

#mu
backward_LRT_mu <- function(fit, alpha = 0.05){

  new_fit <- drop1_LRT_mu(fit, alpha)
  re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit

  #Check if the above new_fit is null model
  if(re_n != 2){
    if(as.character(new_fit$Output_fit$mu)[2] == "1"){
      while (re_n != 2) {
        mu      <- new_fit$Output_fit$mu
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit <- drop1_LRT_mu(fit, alpha)
      }
    }


  }

  return(new_fit)
}


