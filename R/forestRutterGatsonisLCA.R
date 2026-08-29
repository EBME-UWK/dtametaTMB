#' Coupled Forest plot for latent class diagnostic test accuracy meta-analysis
#'
#' Study-specific sensitivities and specificities are empirical Bayes estimates
#' derived from the fitted latent class model. Confidence intervals are based on 
#' the conditional posterior variance of the study-specific random effects 
#' obtained from the TMB Laplace approximation.
#' 
#' @param x Object of class \code{"RutterGatsonisLCA"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest RutterGatsonisLCA
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qnorm plogis qlogis
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.RutterGatsonisLCA <- forest.ReitsmaLCA
