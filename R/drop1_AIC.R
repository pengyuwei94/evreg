
#' Drop one possible covariate on GEV parameter based on AIC
#'
#' Drop a single term to either mu, sigma, and xi based on AIC.
#'
#' drop1_AIC_mu(fit)
#' drop1_AIC_sigma(fit)
#' drop1_AIC_xi(fit)
#'
#' @param fit A model of class "gevreg".
#' @details Add details.
#' @return A list which has the following components
#'     \item{Input_fit}{The input object of the class gevreg.}
#'     \item{Note}{A message that will be printed when input fit and output
#'     fit are the same.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class gevreg if output fit is different
#'     from the input fit.}
#'     \item{AIC}{AIC values for both input model and output model if two models
#'     are different.}
#' @examples
#'
#' ### Annual Maximum and Minimum Temperature
#' # Parameter mu
#' library(extRemes)
#' data(PORTw)
#' PORTw$Year <- (PORTw$Year - min(PORTw$Year)) / (max(PORTw$Year) - min(PORTw$Year))
#' P6 <- gevreg(TMX1, data = PORTw, mu = ~MTMAX + AOindex + Year + STDTMAX + STDMIN + MDTR)
#' drop1_AIC_mu(P6)
#'
#' # Parameter sigma
#' P7 <- gevreg(TMX1, data = PORTw, mu = ~MTMAX + STDTMAX, sigma = ~MTMAX + STDTMAX)
#' drop1_AIC_sigma(P7)
#'
#' # Parameter xi
#'
#'
#' ### Oxford and Worthing annual maximum temperatures
#' #Parameter mu
#' ow$year <- (ow$year - 1901) / (1980 - 1901)
#' ow1 <- gevreg(y = temp, data = ow[-3], mu = ~loc + year, sigma = ~loc,
#' xi = ~loc, sigmalink = identity)
#' drop1_AIC_mu(ow1)
#'
#'
#'
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


    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 in mu
    #    covariate from the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      mu  <- update(fit$formulae$mu, paste("", name[i], sep = "~.-"))  #update mu formula
      new_fit     <- update(fit, mu = mu)     #update fit
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))


    # 3. Identify which of these models has the smallest AIC.
    # 4. If this AIC is smaller than that of the current model then
    #    return this model.  Otherwise, return the original model

    if(min(aic) < AIC(fit)){
      output <- list()
      output$Input_fit <- fit$call

      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$AIC        <- c(AIC(fit), min(aic))
      names(output$AIC) <- c("Input model", "Output model")

      return(output)
    }else{
      output <- list()
      output$Input_fit <- fit$call
      output$Note      <- ("Input fit and output fit are the same.")

      return(output)
    }

  }
  #else for pp fit

}



drop1_AIC_sigma <- function(fit){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if input fit has no covariates on sigma
  if(length(attr(fit$data$D$sigma, "assign")) == 1){
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
    # 1. Fit all of the possible models obtained by dropping 1 covariate on sigma
    #    from the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      sigma  <- update(fit$formulae$sigma, paste("", name[i], sep = "~.-"))  #update sigma formula
      new_fit     <- update(fit, sigma = sigma)     #update fit
      m_list[[i]] <- new_fit
      aic[i]      <- AIC(m_list[[i]])         #get AIC for update fit
    }

    x_i  <- which(aic == min(aic))


    # 3. Identify which of these models has the smallest AIC.
    # 4. If this AIC is smaller than that of the current model then
    #    return this model.  Otherwise, return the original model

    if(min(aic) < AIC(fit)){
      output <- list()
      output$Input_fit <- fit$call

      list <- list()
      list$sigma  <- m_list[[x_i]]$formulae$sigma
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$AIC        <- c(AIC(fit), min(aic))
      names(output$AIC) <- c("Input model", "Output model")

      return(output)
    }else{
      output <- list()
      output$Input_fit <- fit$call
      output$Note      <- ("Input fit and output fit are the same.")

      return(output)
    }

  }
  #else for pp fit

}
