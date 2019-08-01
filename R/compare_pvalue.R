#' Likehood-ratio Test for evreg Objects
#'
#' Compare two evreg objects by using chi-squared distritbuion
#'
#' compare_pvalue(model1, model2)
#'
#' @param model1 A model of class "evreg" which has to be nested within model 2.
#' @param model2 A model of class "evreg".
#' @details Add details.
#' @return A numeric number that is between 0 to 1.
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' f1 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ SOI)
#' compare_pvalue(f0, f1)
#'
#' @export
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
  n1 <- length(name1)
  n2 <- length(name2)
  delta <- 2 * (loglik2 - loglik1)
  p_value <- stats::pchisq(delta, n2 - n1, lower.tail = FALSE)
  return(p_value)
}
