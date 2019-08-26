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

#' Annual Maximum and Minimum Temperature
#'
#' Annual maximum and minimum Winter temperature (degrees centigrade)
#' with a covariate for the North Atlantic Oscillation index from 1927 through 1995.
#' Data is for Winter for Port Jervis, New York (PORTw) and Spring for Sept-Iles, Quebec (SEPTsp).
#'
#' @format A dataframe includes 68 observations and 9 columns
#'   \itemize{
#'     \item{\code{Year}: }{A numeric vector of the calendar year.}
#'     \item{\code{MTMAX}: }{A numeric vector of the means of the winter maximum temperature.}
#'     \item{\code{MTMIN}: }{A numeric vector of the means of the winter minimum temperature.}
#'     \item{\code{STDTMAX}: }{A numeric vector of the standard deviations of the winter maximum temperature.}
#'     \item{\code{STDMIN}: }{A numeric vector of the standard deviations of the winter minimum temperature.}
#'     \item{\code{TMX1}: }{A numeric vector of the annual maximum winter temperature (degrees centigrade).}
#'     \item{\code{TMN0}: }{A numeric vector of the annual minimum winter temperature (degrees centigrade).}
#'     \item{\code{MDTR}: }{A numeric vector of the mean winter diurnal temperature range.}
#'     \item{\code{AOindex}: }{A numeric vector of the North Atlantic Oscillation index (AOindex) from 1927 through 1995}
#'  }
#' @source Eric Gilleland, Richard W. Katz (2016). extRemes 2.0: An Extreme Value Analysis Package
#'     in R. Journal of Statistical Software, 72(8), 1-39. doi:10.18637/jss.v072.i08
"PORTw"

#' Annual Winter Maximum Wave Height
#'
#' Annual Winter Maximum Wave Height in the North Sea from from 1/10/1954 to 31/3/2010,
#' that is, the winters from 56 water years.
#'
#' @format A dataframe contains 56 rows and 6 variables
#'   \itemize{
#'     \item{\code{Hs}: }{A numeric vector of the largest storm maximum value (in metres) observed
#'     over a given winter, which is a measure of sea surface roughness in the northern North Sea.}
#'     \item{\code{waterYear}: }{A numeric vector of the year in which this winter ends
#'     (a water year starts on 1st October)}
#'     \item{\code{meanNAO}: }{A numeric vector of the mean value of the North Atlantic
#'     Oscillation (NAO), which includes climate indice, over a winter.}
#'     \item{\code{maxNAO}: }{A numeric vector of the maximum value of the NAO index over a winter.}
#'     \item{\code{meanAO}: }{A numeric vector of the mean value of the Arctic Oscillation (AO),
#'     which contains climate indice, over a winter.}
#'     \item{\code{maxAO}: }{A numeric vector of the maximum value of the AO index over a winter.}
#'   }
"wm"
