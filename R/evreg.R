#' evreg: Extreme value regression modelling
#'
#' Brief description of the package
#'
#' @details More detail.
#' @seealso \code{\link{gevreg}} GEV generalized linear regression modelling.
#' @docType package
#' @name evreg
#' @import methods
NULL

#' Oxford and Worthing annual maximum temperatures
#'
#' Annual maximum temperatures at Oxford and Worthing (England), for the
#' period 1901 to 1980.
#'
#' @format A dataframe with 80 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{temp}: }{annual maximum temperatures in degrees
#'       Fahrenheit.}
#'     \item{Column 2, \code{year}: }{year in which the maximum was recorded.}
#'     \item{Column 3, \code{name}: }{name of location, "oxford" or "worthing"}
#'     \item{Column 4, \code{loc}: }{location: 1 for "oxford", -1 for
#'       "worthing"}
#'  }
#' @source Tabony, R. C. (1983) Extreme value analysis in meteorology.
#'  \emph{The Meteorological Magazine}, \strong{112}, 77-98.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
"ow"

#' Annual Maximum Sea Levels at Fremantle, Western Australia
#'
#' Annual maximimum sea levels recorded at Fremantle, Western Australia, for
#' the period 1897 to 1989. Also included are annual mean values of the
#' Southern Oscillation Index (SOI), which is a proxy for meteorological
#' volatility in the Southern Hemisphere.
#'
#' @format A dataframe with 86 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{Year}: }{A numeric vector of years..}
#'     \item{Column 2, \code{Sealevel}: }{A numeric vector of annual sea level
#'       maxima..}
#'     \item{Column 3, \code{SOI}: }{A numeric vector of annual mean values of
#'       the Southern Oscillation Index}
#'     \item{Column 4, \code{Year01}: }{A shifted and scaled version of
#'       \code{Year} in which 1987 is mapped to 0 and 1989 to 1.
#'       Scaling of this kind may be useful to avoid numerical difficulties
#'       when using year as a covariate in a regression model.}
#'  }
#' @source Coles, S. G. (2001) An Introduction to Statistical Modelling of
#'   Extreme Values. London: Springer.
"fremantle"
