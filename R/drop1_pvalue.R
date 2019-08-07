#' Drop one possible covariate on GEV parameter based on individual p value
#'
#' Drop a single term to either mu, sigma, and xi based on individual p value.
#'
#' drop1_p_mu(fit, alpha = 0.05)
#' drop1_p_sigma(fit, alpha = 0.05)
#' drop1_p_xi(fit, alpha = 0.05)
#'
#' @param fit A model of class "gevreg".
#' @param alpha Significance level. Default value is 0.05..
#' @details Add details.
#' @return A list which has the following components
#'     \item{Input_fit}{The input object of the class gevreg.}
#'     \item{Note}{A message that will be printed when input fit and output
#'     fit are the same.}
#'     \item{Output_fit}{A list that contains formulae for the parameter,
#'     and the output object of the class gevreg if output fit is different
#'     from the input fit.}
#'     \item{pvalue}{A data frame that contains p value with five decimal
#'     places of the dropped covariate if there is one.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f3 <- gevreg(y = SeaLevel, data = evreg::fremantle[-1], mu = ~Year01 + SOI)
#' drop1_p_mu(f3)
#'
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
  n_mu <- length(attr(fit$data$D$mu, "assign"))
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

    ##likelihood-ratio-test
    p_vec   <- c()
    m_list  <- list()

    #General steps
    # 1. Fit all of the possible models obtained by dropping 1 covariate from the original model
    # 2. Obtain the individual p value of the dropped covariate.

    for(i in 1:length(name)){
      mu          <- update(fit$formulae$mu, paste("", name[i], sep = "~.-")) #update mu formula
      m_list[[i]] <- update(fit, mu = mu)          #update a model call by dropping one covariate on mu

      fit2        <- m_list[[i]]
      p_vec[i]    <- unname(summary(fit2)$pvalue)[n_mu+1] #store all the p-values in one vector
    }
    x_i  <- which(p_vec == max(p_vec))

    ##check significance
    # 3. If none of the p values are significant, return a list that contains three things.
    # Otherwise, return a list of length two.
    if(max(p_vec) > alpha){
      output <- list()
      output$Input_fit <- fit$call

      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$pvalue <- as.data.frame(round(p_vec[x_i],5))
      row.names(output$pvalue) <- name[x_i]
      colnames(output$pvalue)  <- c("Pr(>|z|)")

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


