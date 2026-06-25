#' Forest plot generic
#'
#' @param x Object
#' @param ... Additional arguments
#'
#' @export
forest <- function(x, ...) {
  UseMethod("forest")
}
