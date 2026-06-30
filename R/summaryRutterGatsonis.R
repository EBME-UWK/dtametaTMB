#' Summary Method for RutterGatsonis Objects
#'
#' Provides a concise summary of a fitted \code{"RutterGatsonis"} model.
#'
#' @description
#' This method extracts key components from a fitted HSROC model object
#' returned by \code{\link{fitRutterGatsonis}}. It returns parameter
#' estimates, sensitivity/specificity summaries, and the recovered
#' Reitsma parametrization.
#'
#' @param object An object of class \code{"RutterGatsonis"} as returned by
#'   \code{\link{fitRutterGatsonis}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item \code{estimates} Parameter estimates with standard errors as returned from TMB reported parameters.
#'   \item \code{sensspec} Estimated sensitivity at the specified
#'     specificity, including confidence intervals.
#'   \item \code{Reitsma_recovered} Recovered parameters in the Reitsma parameterization.
#' }
#'
#' @seealso
#' \code{\link{fitRutterGatsonis}}
#'
#' @export
summary.RutterGatsonis <- function(object, ...) {
  return(list(
    estimates = object$sdreport2,
    sensspec = object$sensspec,
    Reitsma_recovered = object$Reitsma_recovered
  ))
}