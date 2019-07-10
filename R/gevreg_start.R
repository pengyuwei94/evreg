gevreg_mustart <- function(data, mulink){
  in2 <- sqrt(6 * stats::var(data$y)) / pi
  in1 <- mean(data$y) - 0.57722 * in2
  return(c(mulink(in1), rep(0, ncol(data$D[[1]]) - 1)))
}

gevreg_sigmastart <- function(data, sigmalink){
  in2 <- sqrt(6 * var(data$y))/pi
  in1 <- mean(data$y) - 0.57722 * in2
  return(c(sigmalink(in2), rep(0, ncol(data$D[[2]]) - 1)))
}

gevreg_xistart <- function(data, xilink){
  return(rep(xilink(0.001), ncol(data$D[[3]])))
}
