#' Summary Method for RutterGatsonisSubgroup Objects
#'
#' Provides a concise summary of a fitted \code{"RutterGatsonisSubgroup"} model.
#'
#' @description
#' This method extracts key components from a fitted HSROC model object
#' returned by \code{\link{fitRutterGatsonisSubgroup}}. It returns parameter
#' estimates and sensitivity/specificity summaries.
#'
#' @param object An object of class \code{"RutterGatsonisSubgroup"} as returned by
#'   \code{\link{fitRutterGatsonisSubgroup}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item \code{estimates} Parameter estimates with standard errors as returned from TMB reported parameters.
#'   \item \code{sensspec} Estimated sensitivity at the specified
#'     specificity, including confidence intervals.
#'   \item \code{Reitsma_recovered} Recovered parameters in the Reitsma parameterization.
#'   \item \code{subgroups} Subgroup names.
#' }
#'
#' @seealso
#' \code{\link{fitRutterGatsonisSubgroup}}
#'
#' @export
summary.RutterGatsonisSubgroup <- function(object, ...) {
  return(list(
    estimates = object$sdreport2,
    sensspec = object$sensspec,
    Reitsma_recovered = object$Reitsma_recovered,
    subgroups = object$subgroups
  ))
}