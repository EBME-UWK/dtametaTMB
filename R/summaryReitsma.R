#' Summary Method for Reitsma Objects
#'
#' Extracts key results from an object of class \code{"Reitsma"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"Reitsma"} as returned by
#'   \code{\link{fitReitsma}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#' }
#'
#' @seealso \code{\link{fitReitsma}}
#'
#'
#' @export
summary.Reitsma <- function(object, ...) {
  return(list(
    estimates = object$estimates,
    sensspec = object$sensspec,
    RutterGatsonis_recovered = object$RutterGatsonis_recovered
  ))
}