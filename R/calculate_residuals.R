gevreg_residuals <- function(x) {
  p <- mu_sigma_xi(x$coefficients, x$data$D)
  shift <- (x$data$y - p[, 1]) / exp(p[, 2])
  return(log1pxdx(shift * p[,3]) * shift)
}

mu_sigma_xi <- function(co, data) {
  np <- length(data)
  p <- vector('list', length = np)
  wh <- 1
  for (i in 1:np){
    which <- wh:(wh - 1 + ncol(data[[i]]))
    p[[i]] <- c(co[which] %*% t(data[[i]]))
    wh <- wh + ncol(data[[i]])
  }
  return(do.call('cbind', p))
}

log1pxdx <- function(x) {
  return(ifelse(x != 0, log1p(x) / x, 1))
}

