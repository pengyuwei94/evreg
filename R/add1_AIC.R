
#' Add one possible covariate on GEV parameter based on AIC
#'
#' Add a single term to either mu, sigma, and xi based on AIC.
#'
#' add1_AIC_mu(fit)
#' add1_AIC_sigma(fit)
#' add1_AIC_xi(fit)
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
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' add1_AIC_mu(f0)
#'
#' @export
add1_AIC_mu <- function(fit){
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
  n_mu <- length(attr(fit$data$D$mu, "assign"))
  if(n_mu == n){
    stop("All covariates has been added on mu for input fit")
  }

  ##identify family from the fit

  ###################################################################
  ####----------------------------GEV----------------------------####
  ###################################################################
  if(inherits(get(m, envir = parent.frame()), "gev")){

    ###Extract covariates on mu if there are in input fit
    if(n_mu != 1){

      x_name <- all.vars(fit$formulae$mu)
      index  <- which(colnames(X) %in% x_name)

      X <- X[-index]   #a data frame with covariates that are not in the sigma formula

    }
    name <- names(X)   #get names from rest of covariates

    #General steps
    # 1. Fit all of the possible models obtained by adding 1 covariate to the original model
    # 2. Calculate the value of AIC for all these models.
    aic <- c()
    m_list <- list()
    for(i in 1:length(name)){
      mu  <- update(fit$formulae$mu, paste("", name[i], sep = "~.+"))  #update mu formula
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
