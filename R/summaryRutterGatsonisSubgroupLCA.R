#' Summary Method for RutterGatsonisSubgroupLCA Objects
#'
#' Extracts key results from an object of class \code{"RutterGatsonisSubgroupLCA"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"RutterGatsonisSubgroupLCA"} as returned by
#'   \code{\link{fitRutterGatsonisSubgroupLCA}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated index test sensitivity at the specified specificity, including confidence intervals.
#'   \item \code{prevref}: Estimated (average) prevalence and reference standard sensitivity/specificitiy with confidence intervals.
#'   \item \code{Reitsma_recovered}: Recovered parameters in the Reitsma parameterization.
#' }
#'
#' @seealso \code{\link{fitRutterGatsonisSubgroupLCA}}
#'
#'
#' @export
summary.RutterGatsonisSubgroupLCA <- function(object, ...) {
  ret <- list(estimates = object$sdreport2,
              sensspec = object$sensspec,
              prevref = object$prevref,
              Reitsma_recovered = object$Reitsma_recovered,
              subgroups = object$subgroups)
  class(ret) <- "summary.RutterGatsonisSubgroupLCA"
  return(ret)
}