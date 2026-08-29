#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#' 
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#'
#' @param x Object of class \code{"RutterGatsonisSubgroup"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param subgroup_label Column name for the subgroup. Defaults to \code{"Subgroup"}.
#' @param order Specifies the ordering of studies
#'   in the forest plot. Can be `"study"` (default), which orders
#'   observations by study identifier and then subgroup, or
#'   `"subgroup"`, which orders observations by subgroup and then
#'   study identifier.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest RutterGatsonisSubgroup
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.RutterGatsonisSubgroup <- forest.ReitsmaSubgroup

