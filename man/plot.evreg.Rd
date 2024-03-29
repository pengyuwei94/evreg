% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_qq_plot.R
\name{plot.evreg}
\alias{plot.evreg}
\title{Plot method for evreg fitted model objects}
\usage{
\method{plot}{evreg}(x, y = NULL, ...)
}
\arguments{
\item{x}{An object returned from either \code{\link{gevreg}} or
\code{\link{ppreg}}.}

\item{y}{Not used.}

\item{...}{Further arguments to be passed to \code{\link[graphics]{plot}}.}
}
\value{
Nothing is returned, only the plot is produced.
}
\description{
Produces a residual QQ plot based on fitted model objects returned from
\code{\link{gevreg}} and \code{\link{ppreg}}.
}
\details{
If \code{x} was returned from \code{\link{gevreg}} then the
 ordered (smallest to largest) residuals \eqn{r_i, i = 1, ..., n} are
 plotted against the 100\eqn{i / (n + 1)}\% \eqn{i = 1, ..., n}
 quantiles of a standard Gumbel distribution.
 If \code{x} was returned from \code{\link{ppreg}} then the
 ordered residuals are plotted against the
 100\eqn{i / (n + 1)}\% \eqn{i = 1, ..., n} of a standard exponential
 distribution.
}
\examples{
# GEV
f1 <- gevreg(SeaLevel, data = evreg::fremantle, mu = ~ SOI)
plot(f1)

# PP (no coavriates)
pp <- ppreg(y = SeaLevel, data = fremantle[, -1])
plot(pp)

# PP regression
pp <- ppreg(y = SeaLevel, data = fremantle[, -1], mu = ~SOI)
}
