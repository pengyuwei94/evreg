#' Add one possible covariate on GEV parameter using individual p value from Wald test
#'
#' Add a single term to either mu, sigma, and xi based on individual p value from Wald test.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05.
#' @details Add details.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if a covariate has been added or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{added_covariate}{A character vector shows added covariate}
#'     \item{pvalue}{A data frame that contains p value with five decimal
#'     places of the added covariate if there is one.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' add1_p_mu(f0)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P0 <- gevreg(y = TMX1, data = PORTw[, -1])
#' add1_p_mu(P0)
#'
#' @name add1_p
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname add1_p
#' @export
add1_p_mu <- function(fit, alpha = 0.05){
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
  n_mu <- ncol(fit$data$D$mu)
  if(n_mu == n){
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
    if(n_mu != 1){
      x_name <- all.vars(fit$formulae$mu)
      index  <- which(colnames(X) %in% x_name)
      X <- X[-index]   #a data frame with covariates that are not in the mu formula
    }
    name <- names(X)   #get names from rest of covariates
    #General steps
    # 1. Fit all of the possible models obtained by adding 1 covariate to the original model
    # 2. Obtain the individual p value of the added covariate.

    p_vec   <- c()
    m_list  <- list()
    for(i in 1:length(name)){
      # adding more variables on mu, one at a time
      mu <- update(fit$formulae$mu, paste("", name[i], sep = "~.+"))  #update mu formula
      # update a model call by adding one additional covariate on mu
      # when we fit a new model in which an extra covariate is added,
      # we use starting values based on the fit of the smaller model.
      # The start value for the new added variable will be set to be zero.
      fcoefs <- fit$coefficients
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- c(unname(fcoefs[1:n_mu]), 0)
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      xistart <- unname(fcoefs[(n_mu + n_sigma + 1):length(fcoefs)])
      m_list[[i]] <- update(fit, mu = mu,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)


      fit2        <- m_list[[i]]
      p_vec[i]  <- unname(summary(fit2)$pvalue)[n_mu+1] #store all the p-values in one vector
    }

    # Check if there is NA values in the p_vec, if so, set it to Infinite
    if(any(is.na(p_vec)) == TRUE){
      na_index <- which(is.na(p_vec) == TRUE)
      p_vec[na_index] <- Inf
    }
    x_i  <- which(p_vec == min(p_vec))

    ##Check significance
    # 3. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if (min(p_vec) < alpha) {
      # Base the output mainly on the chosen fitted model object
      output <- m_list[[x_i]]
      output$added_covariate <- paste0("mu:", name[x_i])
      output$Note <- "covariate added"
      output$Input_fit <- fit$call
      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(p_vec[x_i],5))
      row.names(output$pvalue) <- paste0("mu:", name[x_i])
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    } else {
      output <- fit
      output$added_covariate <- paste0("mu:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
  }
  return(output)
}


# ----------------------------- sigma ---------------------------------

#' @rdname add1_p
#' @export
add1_p_sigma <- function(fit, alpha = 0.05){
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

  ##3. Check if input fit already has all covariates on sigma
  n_sigma <- ncol(fit$data$D$sigma)
  if(n_sigma == n){
    stop("All covariates has been added on sigma for input fit")
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

    ###Extract covariates on sigma if there are in input fit
    if(n_sigma != 1){
      x_name <- all.vars(fit$formulae$sigma)
      index  <- which(colnames(X) %in% x_name)
      X <- X[-index]   #a data frame with covariates that are not in the sigma formula
    }
    name <- names(X)   #get names from rest of covariates
    #General steps
    # 1. Fit all of the possible models obtained by adding 1 covariate to the original model
    # 2. Obtain the individual p value of the added covariate.

    p_vec   <- c()
    m_list  <- list()
    for(i in 1:length(name)){
      # adding more variables on sigma, one at a time
      sigma <- update(fit$formulae$sigma, paste("", name[i], sep = "~.+"))  #update sigma formula
      # update a model call by adding one additional covariate on sigma
      # when we fit a new model in which an extra covariate is added,
      # we use starting values based on the fit of the smaller model.
      # The start value for the new added variable will be set to be zero.
      fcoefs <- fit$coefficients
      n_mu   <- ncol(fit$data$D$mu)
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- c(unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)]), 0)
      xistart <- unname(fcoefs[(n_mu + n_sigma + 1):length(fcoefs)])
      m_list[[i]] <- update(fit, mu = fit$formulae$mu, sigma = sigma,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)


      fit2        <- m_list[[i]]
      p_vec[i]  <- unname(summary(fit2)$pvalue)[n_mu + n_sigma + 1] #store all the p-values in one vector
    }

    # Check if there is NA values in the p_vec, if so, set it to Infinite
    if(any(is.na(p_vec)) == TRUE){
      na_index <- which(is.na(p_vec) == TRUE)
      p_vec[na_index] <- Inf
    }
    x_i  <- which(p_vec == min(p_vec))

    ##Check significance
    # 3. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if (min(p_vec) < alpha) {
      # Base the output mainly on the chosen fitted model object
      output <- m_list[[x_i]]
      output$added_covariate <- paste0("sigma:", name[x_i])
      output$Note <- "covariate added"
      output$Input_fit <- fit$call
      list <- list()
      list$sigma <- m_list[[x_i]]$formulae$sigma
      list$fit   <- m_list[[x_i]]$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(p_vec[x_i],5))
      row.names(output$pvalue) <- paste0("sigma:", name[x_i])
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    } else {
      output <- fit
      output$added_covariate <- paste0("sigma:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
  }
  return(output)
}


# ----------------------------- xi ---------------------------------

#' @rdname add1_p
#' @export
add1_p_xi <- function(fit, alpha = 0.05){
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

  ##3. Check if input fit already has all covariates on xi
  n_xi <- ncol(fit$data$D$xi)
  if(n_xi == n){
    stop("All covariates has been added on xi for input fit")
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

    ###Extract covariates on xi if there are in input fit
    if(n_xi != 1){
      x_name <- all.vars(fit$formulae$xi)
      index  <- which(colnames(X) %in% x_name)
      X <- X[-index]   #a data frame with covariates that are not in the xi formula
    }
    name <- names(X)   #get names from rest of covariates
    #General steps
    # 1. Fit all of the possible models obtained by adding 1 covariate to the original model
    # 2. Obtain the individual p value of the added covariate.

    p_vec   <- c()
    m_list  <- list()
    for(i in 1:length(name)){
      # adding more variables on xi, one at a time
      xi <- update(fit$formulae$xi, paste("", name[i], sep = "~.+"))  #update xi formula
      # update a model call by adding one additional covariate on xi
      # when we fit a new model in which an extra covariate is added,
      # we use starting values based on the fit of the smaller model.
      # The start value for the new added variable will be set to be zero.
      fcoefs  <- fit$coefficients
      n_mu    <- ncol(fit$data$D$mu)
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      xistart <- c(unname(fcoefs[(n_mu + n_sigma + 1):length(fcoefs)]), 0)
      m_list[[i]] <- try(update(fit, mu = fit$formulae$mu, sigma = fit$formulae$sigma,
                                xi = xi,
                               mustart = mustart,
                               sigmastart = sigmastart,
                               xistart = xistart))


      fit2        <- m_list[[i]]
      p_vec[i]  <- unname(summary(fit2)$pvalue)[length(fcoefs)+1] #store all the p-values in one vector
    }

    # Check if there is NA values in the p_vec, if so, set it to Infinite
    if(any(is.na(p_vec)) == TRUE){
      na_index <- which(is.na(p_vec) == TRUE)
      p_vec[na_index] <- Inf
    }
    x_i  <- which(p_vec == min(p_vec))
    ##Check significance
    # 3. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if (min(p_vec) < alpha) {
      # Base the output mainly on the chosen fitted model object
      output <- m_list[[x_i]]
      output$added_covariate <- paste0("xi:", name[x_i])
      output$Note <- "covariate added"
      output$Input_fit <- fit$call
      list <- list()
      list$xi  <- m_list[[x_i]]$formulae$xi
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(p_vec[x_i],5))
      row.names(output$pvalue) <- paste0("xi:", name[x_i])
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    } else {
      output <- fit
      output$added_covariate <- paste0("xi:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
  }
  return(output)
}



