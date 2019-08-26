#' Drop one possible covariate on GEV parameter based on individual p value from Wald test
#'
#' Drop a single term to either mu, sigma, and xi based on individual p value from Wald test.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05.
#' @return An object (a list) of class \code{c("gev", "evreg")} summarising
#'   the new model fit (which may be the same as \code{fit}) and containing the
#'   following additional components
#'     \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
#'     \item{Note}{A message that tells if a covariate has been dropped or not.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class \code{c("gev", "evreg")} if the output fit
#'     is different from the input fit.}
#'     \item{dropped_covariate}{A character vector shows added covariate}
#'     \item{pvalue}{A data frame that contains p value with five decimal
#'     places of the dropped covariate if there is one.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f3 <- gevreg(y = SeaLevel, data = evreg::fremantle[-1], mu = ~Year01 + SOI)
#' drop1_p_mu(f3)
#'
#'
#' ### Annual Maximum and Minimum Temperature
#'
#' P3 <- gevreg(y = TMX1, data = PORTw[, -1], mu = ~MTMAX + STDTMAX + STDMIN)
#' P  <- gevreg(y = TMX1, data = PORTw[, -1], sigma = ~MTMAX + MTMIN + STDTMAX)
#' P5 <- gevreg(y = TMX1, data = PORTw[, -1], xi = ~MTMAX + AOindex)
#' drop1_p_mu(P3)
#' drop1_p_sigma(P)
#' drop1_p_xi(P5)
#' @name drop1_p
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname drop1_p
#' @export
drop1_p_mu <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if there are covariate effects on mu
  n_mu <- ncol(fit$data$D$mu)
  if(n_mu == 1){
    stop("Input fit has no covariate effects on mu")
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
    # General steps
    # 1. drop 1 covariate on mu from the orginal fit according to the p value from Wald test.
    #    Drop the covariate on mu with largest p value if the p value is non-significant.
    p_values  <- summary(fit)$pvalue[2:n_mu]
    # Check if there is NA values in the p_vec, if so, set it to 0
    if(any(is.na(p_values)) == TRUE){
      na_index <- which(is.na(p_values) == TRUE)
      p_values[na_index] <- 0
    }
    index     <- unname(which(p_values == max(p_values)))
    #If there are more then one covariates have max p values,
    #randomly pick one covariate.
    if(length(index) != 1){
      index <- sample(index, 1)
    }
    largest_p <- unname(p_values[index])
    var_name  <- names(p_values[index])
    # check significance
    # 2. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if(largest_p > alpha){
      drop_name <- substring(var_name, first = 5)
      mu <- update(fit$formulae$mu, paste("", drop_name, sep = "~.-"))
      # update a model call by dropping one covariate on mu
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- fcoefs[1:n_mu][-(index+1)]
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma)
      fit2    <- update(fit, mu = mu,
                        sigma = fit$formulae$sigma, xi = fit$formulae$xi,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)
      # Base the output mainly on the chosen fitted model object
      output <- fit2
      output$dropped_covariate <- paste0("mu:", drop_name)
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$mu  <- fit2$formulae$mu
      list$fit <- fit2$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(largest_p, 5))
      row.names(output$pvalue) <- paste0("mu:", drop_name)
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    }else {
      output <- fit
      output$dropped_covariate <- paste0("mu:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
  return(output)
  }
}



# ----------------------------- sigma ---------------------------------

#' @rdname drop1_p
#' @export
drop1_p_sigma <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if there are covariate effects on sigma
  n_sigma <- ncol(fit$data$D$sigma)
  if(n_sigma == 1){
    stop("Input fit has no covariate effects on sigma")
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
    # General steps
    # 1. drop 1 covariate on sigma from the orginal fit according to the p value from Wald test.
    #    Drop the covariate on sigma with largest p value if the p value is non-significant.
    n_mu      <- ncol(fit$data$D$mu)
    p_values  <- summary(fit)$pvalue[(n_mu + 2):(n_mu + n_sigma)]
    # Check if there is NA values in the p_vec, if so, set it to 0
    if(any(is.na(p_values)) == TRUE){
      na_index <- which(is.na(p_values) == TRUE)
      p_values[na_index] <- 0
    }
    index     <- unname(which(p_values == max(p_values)))
    #If there are more then one covariates have max p values,
    #randomly pick one covariate.
    if(length(index) != 1){
      index <- sample(index, 1)
    }
    largest_p <- unname(p_values[index])
    var_name  <- names(p_values[index])
    # check significance
    # 2. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if(largest_p > alpha){
      drop_name <- substring(var_name, first = 8)
      sigma <- update(fit$formulae$sigma, paste("", drop_name, sep = "~.-"))
      # update a model call by dropping one covariate on sigma
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- fcoefs[(n_mu + 1):(n_mu + n_sigma)][-(n_mu + index + 1)]
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma)
      fit2    <- update(fit, sigma = sigma, xi = fit$formulae$xi,
                        mustart = mustart,
                        sigmastart = sigmastart,
                        xistart = xistart)
      # Base the output mainly on the chosen fitted model object
      output <- fit2
      output$dropped_covariate <- paste0("sigma:", drop_name)
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$sigma <- fit2$formulae$sigma
      list$fit   <- fit2$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(largest_p, 5))
      row.names(output$pvalue) <- paste0("sigma:", drop_name)
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    }else {
      output <- fit
      output$dropped_covariate <- paste0("sigma:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
    return(output)
  }
}



# ----------------------------- xi ---------------------------------

#' @rdname drop1_p
#' @export
drop1_p_xi <- function(fit, alpha = 0.05){
  ##1. Check if input arguments are missing
  if(missing(fit))  stop("fit must be specified")

  ##2. Check if input fit is an 'evreg' objects
  m <- deparse(substitute(fit))
  if(!inherits(get(m, envir = parent.frame()), "evreg")){
    stop("Use only with 'evreg' objects")
  }

  ##3. Check if there are covariate effects on xi
  n_xi <- ncol(fit$data$D$xi)
  if(n_xi == 1){
    stop("Input fit has no covariate effects on xi")
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
    # General steps
    # 1. drop 1 covariate on xi from the orginal fit according to the p value from Wald test.
    #    Drop the covariate on xi with largest p value if the p value is non-significant.
    n_mu      <- ncol(fit$data$D$mu)
    n_sigma   <- ncol(fit$data$D$sigma)
    fcoefs    <- fit$coefficients
    p_values  <- summary(fit)$pvalue[(length(fcoefs)-n_xi+1):length(fcoefs)]
    # Check if there is NA values in the p_vec, if so, set it to 0
    if(any(is.na(p_values)) == TRUE){
      na_index <- which(is.na(p_values) == TRUE)
      p_values[na_index] <- 0
    }
    index     <- unname(which(p_values == max(p_values)))
    #If there are more then one covariates have max p values,
    #randomly pick one covariate.
    if(length(index) != 1){
      index <- sample(index, 1)
    }
    largest_p <- unname(p_values[index])
    var_name  <- names(p_values[index])
    # check significance
    # 2. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if(largest_p > alpha){
      drop_name <- substring(var_name, first = 5)
      xi <- update(fit$formulae$xi, paste("", drop_name, sep = "~.-"))
      # update a model call by dropping one covariate on xi
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      mustart <- unname(fcoefs[1:n_mu])
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      # Non-zero components of xistart may mean that the likelihood is zero
      # at the starting values.  This is because for xi not equal to zero
      # there is a constraint on the parameter space.  To avoid this set all
      # components of xistart to 0, i.e. the Gumbel case.
      xistart <- rep(0, length(fcoefs) - n_mu - n_sigma - 1)
      fit2    <- try(update(fit, xi = xi,
                           mustart = mustart,
                           sigmastart = sigmastart,
                           xistart = xistart))
      # Base the output mainly on the chosen fitted model object
      output <- fit2
      output$dropped_covariate <- paste0("xi:", drop_name)
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$xi  <- fit2$formulae$xi
      list$fit <- fit2$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(largest_p, 5))
      row.names(output$pvalue) <- paste0("xi:", drop_name)
      colnames(output$pvalue)  <- c("Pr(>|z|)")
    }else {
      output <- fit
      output$dropped_covariate <- paste0("xi:", "NULL")
      output$Note <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call
    }
    return(output)
  }
}





