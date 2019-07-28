#### mu
#### add 1 covariate on mu based on Likelihood-ratio-test
#### alpha equals 0.05 by default

add1_LRT_mu <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  data    <- eval(fit$call$data)         #save data
  y       <- fit$call$y                  #response variable
  y_index <- which(colnames(data) == y)  #index of column y
  X       <- data[-y_index]              #data only with covariates
  n       <- ncol(X) + 1          #number of covatiates +1 in the data

  ##3. Check if input fit already has all covariates on mu
  if(length(attr(fit$data$D$mu, "assign")) == n){
    stop("All covariates has been added on mu for input fit")
  }

  ##4. Check if input alpha is valid
  if(alpha>1 || alpha <0){
    stop("Significance level alpha should be smaller than 1 and larger than 0")
  }

  ##identify family from the fit

  ###################################################################
  ####----------------------------GEV----------------------------####
  ###################################################################
  if(inherits(get(m, envir = parent.frame()), "gev")){

    ###Extract covariates on mu if there are in input fit
    if(length(attr(fit$data$D$mu, "assign")) != 1){

      x_name <- all.vars(fit$formulae$mu)
      index  <- which(colnames(X) == x_name)

      X <- X[-index]   #a data frame with covariates that are not in the sigma formula

    }
    name <- names(X)   #get names from rest of covariates

    ##likelihood-ratio-test
    p_table <- c()
    m_list  <- list()

    for(i in 1:length(name)){
      mu          <- update(fit$formulae$mu, paste("", name[i], sep = "~.+"))  #update mu formula
      m_list[[i]] <- update(fit, mu = mu)           #update a model call by adding one additional covariate on mu

      fit2        <- m_list[[i]]
      p_table[i]  <- compare_pvalue(fit, fit2)      #store all the p-values in one vector
    }
    x_i  <- which(p_table == min(p_table))

    ##check significance
    if((!(is.nan(p_table))) & (min(p_table) < alpha)){
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





