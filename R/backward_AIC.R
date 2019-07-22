#AIC-based backward elimination

#mu
backward_AIC_mu <- function(fit){

  new_fit <- drop1_AIC_mu(fit, alpha)

  while (as.numeric(length(new_fit)) != 1) {
    mu      <- new_fit$Output_fit$mu
    fit     <- eval(new_fit$Output_fit$fit)
    new_fit <- drop1_AIC_mu(fit)
  }

  return(new_fit)
}
