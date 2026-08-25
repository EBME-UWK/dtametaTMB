#' Forest plot generic
#'
#' Produces coupled forest plots.
#'
#' @param x Object
#' @param ... Additional arguments
#'
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest <- function(x, ...) {
  UseMethod("forest")
}
