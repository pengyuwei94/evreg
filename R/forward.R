#forward

forward <- function(y, data, family){
  #forward selection for gev fit
  if(family == gev){
    x_name <- variable.names(data[-y])
    mu_formula <- as.formula(paste("mu = ~", paste(x_name,collapse = "+")))
    gevreg(y, data, mu_formula)
    #gevreg(SeaLevel, evreg::fremantle, as.formula(paste("mu = ~", paste(variable.names(evreg::fremantle[-c(1,2)]),collapse = "+"))))
    #lm(paste(response,"~",paste(reg_var,collapse = "+")),data)

  }else{# forward selection for pp fit

  }
}
