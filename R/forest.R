#' Forest plot generic
#'
#' Produces coupled forest plots.
#'
#' @param x Object
#' @param ... Additional arguments
#'
#' @return
#' No return value. Called for its side effect of producing a plot.
#' 
#' @export
forest <- function(x, ...) {
  UseMethod("forest")
}
