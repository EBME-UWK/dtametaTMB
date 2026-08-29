#' Coupled Forest plot for latent class diagnostic test accuracy meta-analysis
#'
#' Study-specific sensitivities and specificities are empirical Bayes estimates
#' derived from the fitted latent class model. Confidence intervals are based on 
#' the conditional posterior variance of the study-specific random effects 
#' obtained from the TMB Laplace approximation.
#'
#' @param x Object of class \code{"RutterGatsonisSubgroupLCA"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param subgroup_label Column name for the subgroup. Defaults to \code{"Subgroup"}.
#' @param order Specifies the ordering of studies
#'   in the forest plot. Can be `"study"` (default), which orders
#'   observations by study identifier and then subgroup, or
#'   `"subgroup"`, which orders observations by subgroup and then
#'   study identifier.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest RutterGatsonisSubgroupLCA
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qnorm qlogis plogis
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.RutterGatsonisSubgroupLCA <- forest.ReitsmaSubgroupLCA