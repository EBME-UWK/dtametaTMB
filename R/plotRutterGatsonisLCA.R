#' Plot Results from a Rutter and Gatsonis LCA Model
#'
#' Produces a summary ROC plot for objects of class \code{"RutterGatsonisLCA"}
#' obtained from \code{\link{fitRutterGatsonisLCA}}. Study-specific sensitivities and 
#' specificities are empirical Bayes estimates derived from the fitted latent 
#' class model. 
#'
#' @param x An object of class \code{"RutterGatsonisLCA"}, as returned by
#'   \code{\link{fitRutterGatsonisLCA}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"eb"}{Size proportional to the precision of the empirical Bayes estimates. Default.}
#'    \item{"equal"}{All studies shown with equal size}
#'    \item{"sampsize"}{Size proportional to sample size}
#'  }  
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#'  
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#'   
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown as rectangles.
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates
#'   \item HSROC curve
#' }
#'
#' @references
#' Liu, Y., Chen, Y., & Chu, H. (2015). 
#' A unification of models for meta-analysis of diagnostic accuracy studies without a gold standard. 
#' \emph{Biometrics}, 71(2), 538-547.
#' \doi{10.1111/biom.12264}
#' 
#' Freeman, S. C., Kerby, C. R., Patel, A., Cooper, N. J.,
#' Quinn, T., & Sutton, A. J. (2019).
#' Development of an interactive web-based tool to conduct
#' and interrogate meta-analysis of diagnostic test accuracy studies:
#' MetaDTA.
#' \emph{BMC Medical Research Methodology}, 19, 81.
#' \doi{10.1186/s12874-019-0724-x}
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#' @seealso \code{\link{fitRutterGatsonisLCA}}
#' @method plot RutterGatsonisLCA
#' @export
plot.RutterGatsonisLCA <- function(x, scale=0.02,size=c("eb","equal","sampsize"), 
                                specrange=c(0.7,0.995),
                                main="Diagnostic Test Accuracy Meta-Analysis", ...) {
  size    <- match.arg(size)
  Lambda  <- x$sdreport2["Lambda", "Estimate"]
  beta    <- x$sdreport2["beta","Estimate"]
  roc_points2 <- getROCpoints(Lambda,beta,specrange=specrange)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty = "s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  pct <- getWEIGHTSLCA(x$data,size)
  symbols(x=1-x$data$spec_eb,y=x$data$sens_eb,rectangles=cbind(pct$sp,pct$se)*scale,inches=FALSE,add=TRUE,fg="darkgray")
  #points(x=XP$FPR,y=XP$sens,pch=0,col="darkgray",cex=2)
  # Add the ROC curve
  points(roc_points2, type="l", lwd=2,ann=FALSE)###
  # Add the legend 
  legend("bottomright", 
         bty ="n",
         legend = c(NA,
                    "HSROC curve",
                    "Empirical Bayes estimates"), 
         pch = c(NA,NA,0), 
         lty = c(NA,1,NA), 
         lwd = c(NA,2,NA), 
         col = c(NA,"black","darkgray"))
  invisible(NULL)
}
