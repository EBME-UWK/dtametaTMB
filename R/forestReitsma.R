#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#'
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#' 
#' @param x Object of class \code{"Reitsma"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest Reitsma
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.Reitsma <- function(x,conflevel=0.95, ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  ss <- getForestSensSpec(x=x,conflevel=conflevel)

  XP <- ss$XP
  XP <- XP[order(XP$study), ]
  dt <- XP[,c("study","TP","FP","FN","TN","senslabel","speclabel")]
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  cc <- colnames(dt) 
  colnames(dt) <- c("Study",cc[2:5],ss$senslab,ss$speclab," ",ss$senslab," ",ss$speclab)
  
  getForestPlot(dt=dt,XP=XP)
}

