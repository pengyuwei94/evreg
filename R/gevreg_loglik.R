gevreg_negloglik <- function(pars, data, invmulink, invsigmalink, invxilink) {
  n_pars <- length(pars)
  response_data <- data$y
  mu_mat <- data$D$mu
  sigma_mat <- data$D$sigma
  xi_mat <- data$D$xi
  n_mu <- ncol(mu_mat)
  n_sigma <- ncol(sigma_mat)
  n_xi <- ncol(xi_mat)
  mu_pars <- pars[1:n_mu]
  sigma_pars <- pars[(n_mu + 1):(n_mu + n_sigma)]
  xi_pars <- pars[(n_mu + n_sigma + 1):n_pars]
  mu <- as.vector(invmulink(mu_mat %*% mu_pars))
  sigma <- as.vector(invsigmalink(sigma_mat %*% sigma_pars))
  xi <- as.vector(invxilink(xi_mat %*% xi_pars))
  # Calculate the loglikelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(response_data, loc = mu, scale = sigma,
                           shape = xi, log = TRUE)
  }
  return(-sum(val))
}

