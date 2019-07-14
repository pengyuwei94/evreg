####add 1 covariate on mu based on Likelihood-ratio-test
####alpha equals 0.05 by default

add1_LRT_mu <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  data    <- eval(fit$call$data)         #save data
  y       <- fit$call$y                  #response variable
  y_index <- which(colnames(data) == y)  #index of column y
  X       <- data[-y_index]              #data only with covariates


  m <- deparse(substitute(fit))
  ##identify family from the fit

  ###################################################################
  ####----------------------------GEV----------------------------####
  ###################################################################
  if(inherits(get(a, envir = parent.frame()), "gev")){

    ###Extract covariates on mu if there are in input fit
    if(length(attr(fit$data$D$mu, "assign")) != 1){

      x_name <- all.vars(fit$formulae$mu)
      index  <- which(colnames(data) == x_name)
      #This code gives ERROR since class(x_name) =call
      #Even using eval(x_name), index outputs 'interger(0)'

      X <- X[-index]   #a data frame with covariates that are not in the mu formula

    }
    name <- names(X)   #get names from rest of covariates

    ##likelihood-ratio-test
    p_table <- c()
    m_list  <- list()

    for(i in 1:length(name)){
      mu          <- update(fit$formulae$mu, ~. + name[i]) #update mu formula
      m_list[[i]] <- update(fit, mu = mu)                  #update a model call by adding one additional covariate on mu
      p_table[i]  <- compare_pvalue(fit, add1)         #store all the p-values in one vector
    }

    ##check significance
    if(min(p_table) < alpha){
      x_i  <- which(p_table == min(p_table))
      return(m_list[[1]])
    }else{
      cat("None of the likelihood-ratio-tests are significant")
      return(fit)
    }
  }
  #else for pp fit

}
