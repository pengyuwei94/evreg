% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drop1_gevreg.R
\name{drop1_gevreg}
\alias{drop1_gevreg}
\title{Drop one possible covariate on GEV parameter}
\usage{
drop1_gevreg(fit, parameter = "mu", alpha = 0.05,
  criterion = "pvalue")
}
\arguments{
\item{fit}{An object of class \code{c("gev", "evreg")} returned from
\code{\link{gevreg}} summarising the current model fit.}

\item{parameter}{A specified parameter that needs to drop one covariate.
Equals \code{"mu"} by default.}

\item{alpha}{Significance level if criterion equals pvalue or LRT.
Default value is 0.05.}

\item{criterion}{Either based \code{LRT}(Likelihood ratio test),
\code{AIC}, or \code{pvalue}(by Wald test). Default criterion is \code{pvalue}.}
}
\value{
An object (a list) of class \code{c("gev", "evreg")} summarising
  the new model fit (which may be the same as \code{fit}) and containing the
  following dropitional components
    \item{Input_fit}{The input object of the class \code{c("gev", "evreg")}.}
    \item{Note}{A message that tells if a covariate has been droped or not.}
    \item{Output_fit}{A list that contains formulae for the parameter,
    and the output object of the class \code{c("gev", "evreg")} if the output fit
    is different from the input fit.}
    \item{droped_covariate}{A character vector shows droped covariate}
    \item{criterion_value}{criterion value for if both input model and output model
    are different.}
}
\description{
Drop a single term to either mu, sigma, and xi based on criterion
Likelihood ratio test, AIC, or p value from Wald test
}
\details{
Non-zero components of inital value of xi may mean that the likelihood is zero
at the starting values.  This is because for xi not equal to zero,
there is a constraint on the parameter space.  To avoid this set all
components of initial value of xi to 0, i.e. the Gumbel case.
}
\examples{

}
