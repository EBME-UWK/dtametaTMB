#' Summary Method for RutterGatsonisReg Objects
#'
#' Provides a concise summary of a fitted \code{"RutterGatsonisReg"} model.
#'
#' @description
#' This method extracts key components from a fitted HSROC model object
#' returned by \code{\link{fitRutterGatsonisReg}}. It returns parameter
#' estimates, sensitivity/specificity summaries.
#'
#' @param object An object of class \code{"RutterGatsonisReg"} as returned by
#'   \code{\link{fitRutterGatsonisReg}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item \code{estimates} Parameter estimates with standard errors as returned from TMB reported parameters.
#'   \item \code{sensspec} Estimated sensitivity at the specified
#'     specificity, including confidence intervals.
#' }
#'
#' @seealso
#' \code{\link{fitRutterGatsonisReg}}
#'
#' @export
summary.RutterGatsonisReg <- function(object, ...) {
  return(list(
    estimates = object$sdreport2,
    sensspec = object$sensspec
  ))
}