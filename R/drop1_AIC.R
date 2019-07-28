#### mu
#### drop 1 covariate on mu based on AIC

drop1_AIC_mu <- function(fit){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if input fit is null model
  if(length(attr(fit$data$D$mu, "assign")) == 1){
    stop("Input fit has no covariate effects on mu")
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

    ###Extract covariates on mu
    x_name <- all.vars(fit$formulae$mu)
    index  <- which(colnames(X) %in% x_name)

    X <- X[index]   #a data frame with covariates that are in the mu formula
    name <- names(X)

    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      mu  <- update(fit$formulae$mu, paste("", name[i], sep = "~.-"))  #update mu formula
      new_fit     <- update(fit, mu = mu)     #update fit
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))

    ##return the fit with smallest AIC
    ##if AIC of dropping one covariate on mu is smaller than AIC of the input fit
    if(min(aic) < AIC(fit)){
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
