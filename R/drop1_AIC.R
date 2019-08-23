#' Drop one possible covariate on GEV parameter based on AIC
#'
#' Drop a single term to either mu, sigma, and xi based on AIC.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if a covariate has been added or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{dropped_covariate}{A character vector shows dropped covariate}
#'     \item{AIC}{AIC values for both input model and output model if two models
#'     are different.}
#' @examples
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P6 <- gevreg(TMX1, data = PORTw[,-1], mu = ~MTMAX + AOindex + STDTMAX + STDMIN + MDTR)
#' P7 <- gevreg(TMX1, data = PORTw[,-1], sigma = ~MTMAX + STDTMAX + STDMIN + MDTR)
#' P8 <- gevreg(TMX1, data = PORTw[,-1], xi = ~MTMAX + STDTMAX + STDMIN + MDTR)
#' drop1_AIC_mu(P6)
#' drop1_AIC_sigma(P7)
#' drop1_AIC_xi(P8)
#'
#'
#' ### Oxford and Worthing annual maximum temperatures
#' #Parameter mu
#' ow$year <- (ow$year - 1901) / (1980 - 1901)
#' ow1 <- gevreg(y = temp, data = ow[-3], mu = ~loc + year, sigma = ~loc,
#' xi = ~loc, sigmalink = identity)
#' drop1_AIC_mu(ow1)
#'
#' @name drop1_AIC
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname drop1_AIC
#' @export
drop1_AIC_mu <- function(fit){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if input fit has no covariates on mu
  n_mu <- ncol(fit$data$D$mu)
  if(n_mu == 1){
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

    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 covariate from the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      #dropping more variables on mu, one at a time
      mu <- update(fit$formulae$mu, paste("", name[i], sep = "~.-")) #update mu formula
      # update a model call by dropping one covariate on mu
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      which_to_drop <- paste0("mu: ", name[i])
      drop_num <- which(which_to_drop == names(fcoefs))
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- fcoefs[1:n_mu][-drop_num]
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma)
      new_fit     <- update(fit, mu = mu,
                            sigma = fit$formulae$sigma, xi = fit$formulae$xi,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))


    # 3. Identify which of these models has the smallest AIC.
    # 4. If this AIC is smaller than that of the current model then
    #    return this model.  Otherwise, return the original model

    if(min(aic) < AIC(fit)){
      output <- m_list[[x_i]]
      output$dropped_covariate <- paste0("mu:", name[x_i])
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$AIC        <- c(AIC(fit), min(aic))
      names(output$AIC) <- c("Input model", "Output model")

    }else{
      output <- fit
      output$dropped_covariate <- paste0("mu:", "NULL")
      output$Note      <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }

  }
  return(output)

}



# ----------------------------- sigma ---------------------------------

#' @rdname drop1_AIC
#' @export
drop1_AIC_sigma <- function(fit){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if input fit has no covariates on sigma
  n_sigma <- ncol(fit$data$D$sigma)
  if(n_sigma == 1){
    stop("Input fit has no covariate effects on sigma")
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

    ###Extract covariates on sigma
    x_name <- all.vars(fit$formulae$sigma)
    index  <- which(colnames(X) %in% x_name)

    X <- X[index]   #a data frame with covariates that are in the sigma formula
    name <- names(X)

    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 covariate from the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      #dropping more variables on sigma, one at a time
      sigma <- update(fit$formulae$sigma, paste("", name[i], sep = "~.-")) #update sigma formula
      # update a model call by dropping one covariate on sigma
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      which_to_drop <- paste0("sigma: ", name[i])
      drop_num <- which(which_to_drop == names(fcoefs))
      n_mu <- ncol(fit$data$D$mu)
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- fcoefs[(n_mu + 1):(n_mu + n_sigma)][-drop_num]
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma)
      new_fit     <- update(fit, sigma = sigma, xi = fit$formulae$xi,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))


    # 3. Identify which of these models has the smallest AIC.
    # 4. If this AIC is smaller than that of the current model then
    #    return this model.  Otherwise, return the original model

    if(min(aic) < AIC(fit)){
      output <- m_list[[x_i]]
      output$dropped_covariate <- paste0("sigma:", name[x_i])
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$sigma <- m_list[[x_i]]$formulae$sigma
      list$fit   <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$AIC        <- c(AIC(fit), min(aic))
      names(output$AIC) <- c("Input model", "Output model")

    }else{
      output <- fit
      output$dropped_covariate <- paste0("sigma:", "NULL")
      output$Note      <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }

  }
  return(output)

}


# ----------------------------- xi ---------------------------------

#' @rdname drop1_AIC
#' @export
drop1_AIC_xi <- function(fit){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if input fit has no covariates on xi
  n_xi <- ncol(fit$data$D$xi)
  if(n_xi == 1){
    stop("Input fit has no covariate effects on xi")
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

    ###Extract covariates on xi
    x_name <- all.vars(fit$formulae$xi)
    index  <- which(colnames(X) %in% x_name)

    X <- X[index]   #a data frame with covariates that are in the xi formula
    name <- names(X)

    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 covariate from the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      #dropping more variables on xi, one at a time
      xi <- update(fit$formulae$xi, paste("", name[i], sep = "~.-")) #update xi formula
      # update a model call by dropping one covariate on xi
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      n_mu <- ncol(fit$data$D$mu)
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma - 1)
      new_fit     <- try(update(fit, xi = xi,
                               mustart = mustart,
                               sigmastart = sigmastart,
                               xistart = xistart))
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))


    # 3. Identify which of these models has the smallest AIC.
    # 4. If this AIC is smaller than that of the current model then
    #    return this model.  Otherwise, return the original model

    if(min(aic) < AIC(fit)){
      output <- m_list[[x_i]]
      output$dropped_covariate <- paste0("xi:", name[x_i])
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$xi  <- m_list[[x_i]]$formulae$xi
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$AIC        <- c(AIC(fit), min(aic))
      names(output$AIC) <- c("Input model", "Output model")

    }else{
      output <- fit
      output$dropped_covariate <- paste0("xi:", "NULL")
      output$Note      <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }

  }
  return(output)

}

