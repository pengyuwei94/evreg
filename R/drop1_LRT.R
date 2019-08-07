#' Drop one possible covariate on GEV parameter based on Likelihood-ratio-test
#'
#' Drop a single term to either mu, sigma, and xi based on likelihood-ratio-test.
#'
#' @param fit An object of class \code{c("gev", "evreg")} returned from
#'   \code{\link{gevreg}} summarising the current model fit.
#' @param alpha Significance level. Default value is 0.05..
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
#'     \item{pvalue}{A data frame that contains p value with five decimal
#'     places of the Likelihood-ratio-test.}
#' @examples
#'
#' ### Oxford and Worthing annual maximum temperatures
#'
#' ow$year <- (ow$year - 1901) / (1980 - 1901)
#' ow1 <- gevreg(y = temp, data = ow[-3], mu = ~loc + year, sigma = ~loc,
#' xi = ~loc, sigmalink = identity)
#' drop1_LRT_mu(ow1)
#' @name drop1_LRT
NULL
## NULL

# ----------------------------- mu ---------------------------------

#' @rdname drop1_LRT
#' @export
drop1_LRT_mu <- function(fit, alpha = 0.05){
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
    index  <- which(colnames(X) %in% x_name)

    X <- X[index]   #a data frame with covariates that are in the mu formula
    name <- names(X)

    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 covariate from the original model
    # 2. Calculate the p value between new model and original model.

    ##likelihood-ratio-test
    p_vec   <- c()
    m_list  <- list()
    for(i in 1:length(name)){
      #dropping more variables on mu, one at a time
      mu          <- update(fit$formulae$mu, paste("", name[i], sep = "~.-")) #update mu formula
      # update a model call by dropping one covariate on mu
      # when we fit a new model in which an covariate is dropped,
      # we use starting values based on the fit of the larger model.
      fcoefs  <- fit$coefficients
      n_sigma <- ncol(fit$data$D$sigma)
      mustart <- c((unname(fcoefs[1:i])), unname(fcoefs[((i+2):n_mu)]))
      sigmastart <- unname(fcoefs[(n_mu + 1):(n_mu + n_sigma)])
      xistart <- unname(fcoefs[(n_mu + n_sigma + 1):length(fcoefs)])
      m_list[[i]] <- update(fit, mu = mu,
                            mustart = mustart,
                            sigmastart = sigmastart,
                            xistart = xistart)


      fit2        <- m_list[[i]]
      p_vec[i]    <- compare_pvalue(fit2, fit)      #store all the p-values in one vector
    }
    x_i  <- which(p_vec == max(p_vec))

    ##check significance
    # 3. If none of the p values are significant, return a new fitted object.
    #    Otherwise, return the old fitted object.
    if(max(p_vec) > alpha){
      output <- m_list[[x_i]]
      output$dropped_covariate <- name[x_i]
      output$Note <- "covariate dropped"
      output$Input_fit <- fit$call
      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list
      output$pvalue <- as.data.frame(round(p_vec[x_i],5))
      row.names(output$pvalue) <- name[x_i]
      colnames(output$pvalue)  <- c("LRT_pvalue")

    }else{
      output <- fit
      output$dropped_covariate <- NULL
      output$Note      <- "Input fit and output fit are the same"
      output$Input_fit <- fit$call

    }
  }
  return(output)
}

# ----------------------------- sigma ---------------------------------


# ----------------------------- xi ---------------------------------


