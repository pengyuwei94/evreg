#AIC-based backward elimination

#mu
backward_AIC_mu <- function(fit){

  new_fit <- drop1_AIC_mu(fit)
  re_n    <- as.numeric(length(new_fit))   #length of return list from new_fit

  #Check if the above new_fit is null model
  if(re_n != 2){
    if(as.character(new_fit$Output_fit$mu)[2] == "1"){
      while (re_n != 2) {
        mu      <- new_fit$Output_fit$mu
        fit     <- eval(new_fit$Output_fit$fit)
        new_fit <- drop1_AIC_mu(fit)
      }
    }


  }

  return(new_fit)
}
