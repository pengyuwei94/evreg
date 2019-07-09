compare_pvalue <- function(model1, model2) {

  #Check that models are specified
  if(missing(model1)) stop("model one must be specified")
  if(missing(model2)) stop("model two must be specified")

  #Check that models are 'evreg' objects
  m1 <- deparse(substitute(model1))
  m2 <- deparse(substitute(model2))
  m  <- c(m1, m2)
  n  <- length(m)
  for(i in 1:n){
    if(!inherits(get(m[i], envir = parent.frame()), "evreg"))
    stop("Use only with 'evreg' objects")
  }

  #Check that model 2 has more parameters than model 1
  if((!all(names(model2$coefficients) %in% names(model1$coefficients)))){
      warning("models may not be nested")
  }

  #Check that the maximised log-likelihoods are the correct way round

  p_value <- stats::pchisq(model1$loglik, model2$loglik)
  return(p_value)
}
