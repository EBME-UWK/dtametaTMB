#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#'
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#' 
#' @param x Object of class \code{"RutterGatsonis"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest RutterGatsonis
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.RutterGatsonis <- forest.Reitsma
