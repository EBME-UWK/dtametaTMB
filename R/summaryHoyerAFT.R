#' Summary Method for HoyerAFT Objects
#'
#' Extracts and returns the summary information stored in a fitted
#' \code{HoyerAFT} model object.
#'
#' This method is part of the standard S3 \code{summary()} generic and
#' provides access to summary statistics computed from the fitted model,
#' as stored in the \code{sdreport2} and \code{sensspec} component of the object.
#'
#' @param object An object of class \code{"HoyerAFT"}, typically the result
#'   of a call to \code{\link{fitHoyerAFT}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with two data frames containing summary statistics derived
#'   from the fitted model, corresponding to the \code{sdreport2} and \code{sensspec} element of the object.
#'
#' @seealso \code{\link{fitHoyerAFT}}, \code{\link[=summary.HoyerAFT]{summary}}
#'
#' @export
summary.HoyerAFT <- function(object, ...) {
  return(list(sdreport2=object$sdreport2,
              sensspec=object$sensspec))
}