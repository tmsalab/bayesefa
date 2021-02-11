#' International Exchange Rate Dataset
#'
#' This dataset was considered by West and Harrison (1997) and 
#' Lopes and West (2004). The dataset consists of n = 143 monthly 
#' first-differences of the exchange rates of 6 international 
#' currencies against the British pound, from Jan 1975 to Dec 1986, 
#' these currencies are: US dollar (USD), Canadian dollar (CAD), 
#' Japanese yen (JPY), French franc (FRF), German (deutsche) mark (DEM),
#' and the Italian lira (ITL).
#'
#' @format A 143 by 7 \code{matrix} of exchange rate time-series. The variables include: 
#' \describe{
#'  \item{\code{Month_Year}}{Month and year of exchange rate data.}
#'  \item{\code{USD}}{US dollar}
#'  \item{\code{CAD}}{Canadian dollar}
#'  \item{\code{JPY}}{Japanese yen}
#'  \item{\code{FRF}}{French franc}
#'  \item{\code{DEM}}{German (deutsche) mark}
#'  \item{\code{ITL}}{Italian lira}
#' }
#' @source 
#' 
#' Lopes, H. F., and West, M. (2004). Bayesian model assessment in factor analysis, Statistica Sinica, 14, 41–67.
#' 
#' Man, A. & Culpepper, S. A. (2020). A mode-jumping algorithm for Bayesian factor analysis. Journal of the American Statistical Association, doi:10.1080/01621459.2020.1773833.
#' 
#' West, M., and Harrison, J. (1997), Bayesian forecasting and dynamic models (2nd ed.), Berlin, Heidelberg: Springer-Verlag.
#' 
#' @author Steven Culpepper
#' 
"exchangerate"
