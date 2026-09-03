#' Summary Method for ReitsmaSubgroupLCA Objects
#'
#' Extracts key results from an object of class \code{"ReitsmaSubgroupLCA"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"ReitsmaSubgroupLCA"} as returned by
#'   \code{\link{fitReitsmaSubgroupLCA}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#' }
#'
#' @seealso \code{\link{fitReitsmaSubgroupLCA}}
#'
#'
#' @export
summary.ReitsmaSubgroupLCA <- function(object, ...) {
  ret <- list(estimates = object$sdreport2,
              sensspec = object$sensspec,
              RutterGatsonis_recovered = object$RutterGatsonis_recovered,
              subgroups = object$subgroups)
  class(ret) <- "summary.ReitsmaSubgroupLCA"
  return(ret)
}