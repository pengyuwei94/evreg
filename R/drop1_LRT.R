#### mu
#### drop 1 covariate on mu based on Likelihood-ratio-test
#### alpha equals 0.05 by default

drop1_LRT_mu <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if there are covariate effects on mu
  if(length(attr(fit$data$D$mu, "assign")) == 1){
    stop("Input fit has no covariate effects on mu")
  }

  ##4. Check if input alpha is valid
  if(alpha>1 || alpha <0){
    stop("Significance level alpha should be smaller than 1 and larger than 0")
  }

  data    <- eval(fit$call$data)         #save data
  y       <- fit$call$y                  #response variable
  y_index <- which(colnames(data) == y)  #index of column y
  X       <- data[-y_index]              #data only with covariates

  ##identify family from the fit

  ###################################################################
  ####----------------------------GEV----------------------------####
  ###################################################################
  if(inherits(get(m, envir = parent.frame()), "gev")){

    ###Extract existing covariates on mu
    x_name <- all.vars(fit$formulae$mu)
    index  <- which(colnames(X) == x_name)

    X <- X[index]   #a data frame with covariates that are in the mu formula
    name <- names(X)

    ##likelihood-ratio-test
    p_table <- c()
    m_list  <- list()

    for(i in 1:length(name)){
      mu          <- update(fit$formulae$mu, paste("", name[i], sep = "~.-")) #update mu formula
      m_list[[i]] <- update(fit, mu = mu)           #update a model call by dropping one covariate on mu

      fit2        <- m_list[[i]]
      p_table[i]  <- compare_pvalue(fit2, fit)      #store all the p-values in one vector
    }
    x_i  <- which(p_table == min(p_table))

    ##check significance
    if(min(p_table) < alpha){
      output <- list()
      output$Input_fit <- fit$call

      list <- list()
      list$mu <- mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      return(output)

    }else{
      output <- list()
      output$Input_fit <- fit$call
      print("Input fit and output fit are the same.")

      return(output)
    }
  }
  #else for pp fit

}



