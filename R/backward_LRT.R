#Likelihood-ratio-test-based backward elimination

#mu
backward_LRT_mu <- function(fit, alpha = 0.05){

  new_fit = drop1_LRT_mu(fit, alpha)

  while (as.numeric(length(new_fit)) != 1) {
    mu      <- new_fit$Output_fit$mu
    fit     <- eval(new_fit$Output_fit$fit)
    new_fit <- drop1_LRT_mu(fit, alpha)
  }

  return(new_fit)
}


