#' Coupled Forest plot for latent class diagnostic test accuracy meta-analysis
#'
#' Study-specific sensitivities and specificities are empirical Bayes estimates
#' derived from the fitted latent class model. Confidence intervals are based on 
#' the conditional posterior variance of the study-specific random effects 
#' obtained from the TMB Laplace approximation.
#' 
#' @param x Object of class \code{"ReitsmaLCA"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest ReitsmaLCA
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qnorm plogis qlogis
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.ReitsmaLCA <- function(x,conflevel=0.95, ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  ss <- getForestSensSpecLCA(x=x,conflevel=conflevel)
  
  XP <- ss$XP
  XP <- XP[order(XP$study), ]
  dt <- XP[,c("study","y11","y10","y01","y00","senslabel","speclabel")]
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  cc <- colnames(dt) 
  colnames(dt) <- c("Study","+/+","+/-","-/+","-/-",ss$senslab,ss$speclab," ",ss$senslab," ",ss$speclab)
  
  getForestPlot(dt=dt,XP=XP)
}