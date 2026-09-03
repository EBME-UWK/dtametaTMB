#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#'
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#' 
#' @param x Object of class \code{"HoyerAFT"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest HoyerAFT
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.HoyerAFT <- function(x,conflevel=0.95, ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  ss <- getForestSensSpec(x=x,conflevel=conflevel)
  
  XP <- ss$XP
  XP <- XP[order(XP$study, XP$threshold), ]
  dt <- XP[, c("study", "threshold", "TP","FP","FN","TN","senslabel", "speclabel")]
  dt$study <- as.character(dt$study)
  dt$study[duplicated(dt$study)] <- " "
    
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  cc <- colnames(dt) 
  colnames(dt) <- c("Study","Threshold",cc[3:6],
                    ss$senslab,ss$speclab," ",ss$senslab," ",ss$speclab)
  
  getForestPlotSub(dt=dt,XP=XP,s=2)
}

