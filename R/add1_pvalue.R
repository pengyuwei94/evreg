#' Add one possible covariate on GEV parameter using individual p value
#'
#' Add a single term to either mu, sigma, and xi based on individual p value.
#'
#' add1_p_mu(fit, alpha = 0.05)
#' add1_p_sigma(fit, alpha = 0.05)
#' add1_p_xi(fit, alpha = 0.05)
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
#'     \item{pvalue}{A p value with five decimal places of the covariate
#'     that is been added if there is one.}
#' @examples
#'
#' ### Fremantle sea levels
#'
#' f0 <- gevreg(SeaLevel, data = evreg::fremantle)
#' add1_p_mu(f0)
#'
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
  n_mu <- length(attr(fit$data$D$mu, "assign"))
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

      X <- X[-index]   #a data frame with covariates that are not in the sigma formula

    }
    name <- names(X)   #get names from rest of covariates

    #General steps
    # 1. Fit all of the possible models obtained by adding 1 covariate to the original model
    # 2. Obtain the individual p value of the added covariate.

    ##Likelihood-ratio-test
    p_table <- c()
    m_list  <- list()

    # Preparation for updating initial values.
    # Detail explanations are in the for loop.
    #n_sigma <- length(attr(fit$data$D$sigma, "assign"))
    #n_coeff <- length(fit$coefficients)

    #muint   <- c(unname(fit$coefficients[1:n_mu]), 0)
    #sigint  <- unname(fit$coefficients[(n_mu + 1):(n_mu + n_sigma)])
    #xiint   <- unname(fit$coefficients[(n_mu + n_sigma + 1):n_coeff])

    #muint   <- deparse(substitute(muint))
    #sigint  <- deparse(substitute(sigint))
    #xiint   <- deparse(substitute(xiint))



    for(i in 1:length(name)){
      # adding more variables on mu, one at a time
      mu          <- update(fit$formulae$mu, paste("", name[i], sep = "~.+"))  #update mu formula

      # update a model call by adding one additional covariate on mu
      # when we fit a new model in which an extra covariate is added,
      # we use starting values based on the fit of the smaller model.
      # The start value for the new added variable will be set to be zero.
      m_list[[i]] <- update(fit, mu = mu
                            #,
                            #mustart = muint,
                            #sigmastart = sigint,
                            #xistart = xiint
      )


      fit2        <- m_list[[i]]
      p_table[i]  <- unname(summary(fit2)$pvalue)[n_mu+1] #store all the p-values in one vector
    }
    x_i  <- which(p_table == min(p_table))

    ##Check significance
    # 3. If none of the p values are significant, return a list that contains three things.
    # Otherwise, return a list of length two.
    if(min(p_table) < alpha){
      output <- list()
      output$Input_fit <- fit$call

      list <- list()
      list$mu  <- m_list[[x_i]]$formulae$mu
      list$fit <- m_list[[x_i]]$call
      output$Output_fit <- list

      output$pvalue <- round(p_table[x_i],5)

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





