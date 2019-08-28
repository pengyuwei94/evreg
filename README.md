
<!-- README.md is generated from README.Rmd. Please edit that file -->
evreg
=====

Extreme Value Regression
------------------------

### What does evreg do?

The `evreg` package provides functions that model non-stationary extreme values using univariate extreme value regression modelling, while also performing variable selection to objectively determine the best model for a given dataset. This is achieved by using generalized linear modelling for each parameter, and by fitting models with maximum likelihood estimates.

### A simple example

One of the main functions in the evreg package is `gevreg`, which performs the model fitting for generalized extreme value regression modelling. Then, `forward_gevreg` can be used to perform forward selection to objectively determine the best model. The following code fits a GEV regression model using the `gevreg` function for `fremantle` dataset and then uses `forward_gevreg` perform forward selection.

``` r
library(evreg)
data <- fremantle[ , -which(names(fremantle) %in% c("Year"))]
fit0 <- gevreg(SeaLevel, data)
forward_gevreg(fit0)
#> 
#> Call:
#> gevreg(y = SeaLevel, data = data, mu = mu, mustart = mustart, 
#>     sigmastart = sigmastart, xistart = xistart)
#> 
#> Convergence: TRUE 
#> 
#> Coefficients:
#>    mu: (Intercept)          mu: Year01             mu: SOI  
#>            1.38433             0.19449             0.05452  
#> sigma: (Intercept)     xi: (Intercept)  
#>           -2.11418            -0.14999  
#> 
#> Log-likelihood: 53.9     DF: 5   AIC: -97.8
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("evreg")
```

### Vignette
