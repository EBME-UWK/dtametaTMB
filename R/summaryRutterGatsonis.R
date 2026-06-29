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
#'   \item \code{estimates} Matrix of reported parameter estimates and standard errors
#'     from the TMB \code{sdreport} or summary (i.e., \code{sdreport2}).
#'   \item \code{sensspec} Data frame with estimated sensitivity at the specified
#'     specificity, including confidence intervals.
#'   \item \code{Reitsma_recovered} Data frame of parameters transformed into the
#'     Reitsma (bivariate) parametrization.
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