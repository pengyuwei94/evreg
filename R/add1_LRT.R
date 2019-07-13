####add 1 covariate on mu based on Likelihood-ratio-test
####alpha equals 0.05 by default

add1_LRT_mu <- function(fit, data, alpha = 0.05){
  #1. Check if fit is missing
  if(missing(fit)) stop("fit must be specified")

  y       <- fit$call$y
  y_index <- which(colnames(data) == y)
  X       <- data[-y_index]


  m <- deparse(substitute(fit))
  ##identify family from the fit

  ###################################################################
  ####----------------------------GEV----------------------------####
  ###################################################################
  if(inherits(get(a, envir = parent.frame()), "gev")){

    ###Extract covariates on mu if there are in input fit
    if(length(attr(fit$data$D$mu, "assign")) != 1){

      x_name <- all.vars(fit$formulae$mu)
      index  <- which(colnames(data) == x_name) #This code gives ERROR since class(x_name) =call
      #a data frame with covariates that are not in the mu formula
      X <- X[-index]

    }
    #get names from rest of covariates
    name <- names(X)
    mu <- fit$formulae$mu
    sigma <- fit$formulae$sigma
    xi <- fit$formulae$xi

    #likelihood-ratio-test
    p_table <- c()
    for(i in 1:length(name)){
      add1 <- gevreg(y, data, mu = as.formula(paste(mu, paste(name[i], collapse = "+"))), sigma = sigma, xi = xi)
      #store all the p-values in one vector
      p_table[i] <- compare_pvalue(fit, add1)
    }

    #check significance
    if(min(p_table) < alpha){
      x_i <- which(p_table == min(p_table))
      add1 <- gevreg(y, data, mu = as.formula(paste(mu, paste(name[xi], collapse = "+"))), sigma = sigma, xi = xi)
      return(add1)
    }else{
      cat("None of the likelihood-ratio-tests are significant")
      return(fit)
    }
  }
  #else for pp fit

}
