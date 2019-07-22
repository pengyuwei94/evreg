#Likehood-ratio test
#Compare two evreg objects by using chi-squared distritbuion
#Model 1 has to be nested within model 2

compare_pvalue <- function(model1, model2) {

  ##1. Check that models are specified
  if(missing(model1)) stop("model one must be specified")
  if(missing(model2)) stop("model two must be specified")

  ##2. Check that models are 'evreg' objects
  m1 <- deparse(substitute(model1))
  m2 <- deparse(substitute(model2))
  m  <- c(m1, m2)
  n  <- length(m)
  for(i in 1:n){
    if(!inherits(get(m[i], envir = parent.frame()), "evreg"))
    stop("Use only with 'evreg' objects")
  }

  ##3. Check that model 2 has more parameters than model 1
  name1 <- names(model1$coefficients)
  name2 <- names(model2$coefficients)
  if((!all(name1 %in% name2))){
      warning("make sure model1 is nested within model2")
  }

  ##4. Check that the maximised log-likelihoods are the correct way round
  loglik1 <- model1$loglik
  loglik2 <- model2$loglik
  if(loglik2 < loglik1) stop("models may not be nested")

  ##Compute p-value
  p_value <- stats::pchisq(loglik1, loglik2)
  return(p_value)
}
