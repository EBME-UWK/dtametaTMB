#' Summary Method for ReitsmaSubgroup Objects
#'
#' Extracts key results from an object of class \code{"ReitsmaSubgroup"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"ReitsmaSubgroup"} as returned by
#'   \code{\link{fitReitsmaSubgroup}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{subgroups} Subgroup names.
#' }
#'
#' @seealso \code{\link{fitReitsmaSubgroup}}
#'
#'
#' @export
summary.ReitsmaSubgroup <- function(object, ...) {
  rnu <- grep("^nu_", rownames(object$estimates_nu))
  nu <- object$estimates_nu[rnu,]
  return(list(
    estimates = rbind(object$estimates_mu,nu),
    sensspec = object$sensspec,
    RutterGatsonis_recovered = object$RutterGatsonis_recovered,
    subgroups = object$subgroups
  ))
}