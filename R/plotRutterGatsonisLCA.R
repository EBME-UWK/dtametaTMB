#' Plot Results from a Rutter and Gatsonis LCA Model
#'
#' Produces a summary ROC plot for objects of class \code{"RutterGatsonisLCA"}
#' obtained from \code{\link{fitRutterGatsonisLCA}}. Study-specific sensitivities and 
#' specificities are empirical Bayes estimates derived from the fitted latent 
#' class model. 
#'
#' @param x An object of class \code{"RutterGatsonisLCA"}, as returned by
#'   \code{\link{fitRutterGatsonisLCA}}.
#' @param symbol Character string specifying the symbol used for study-level
#'   estimates. One of \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. The default is \code{"rectangle"}.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally.
#' @param size Character string controlling display of study-level
#'   point estimates:
#'   \describe{
#'     \item{\code{"eb"}}{
#'       Horizontal and vertical symbol dimensions are proportional to
#'       normalized inverse posterior standard errors of the
#'       study-specific empirical Bayes logit-specificity and
#'       logit-sensitivity estimates, respectively.
#'     }
#'     \item{\code{"equal"}}{
#'       All studies are shown with equal size.
#'     }
#'     \item{\code{"sampsize"}}{
#'       Horizontal and vertical symbol dimensions are both proportional
#'       to relative total study sample size. The same dimension measure
#'       is used horizontally and vertically because true disease status
#'       is latent.
#'     }
#'     \item{\code{"eb_revman"}}{
#'       RevMan-style graphical scaling applied to the precision of the
#'       empirical Bayes estimates. This option uses the RevMan-style
#'       fixed baseline and point-scaling constants but is not an
#'       original RevMan latent-class weighting method.
#'     }
#'     \item{\code{"sampsize_revman"}}{
#'       RevMan-style graphical scaling based on total study sample size.
#'       The same transformed sample-size measure is used for the
#'       horizontal and vertical symbol dimensions.
#'     }
#'   }
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
#' Study-specific estimates are shown using the selected symbol (rectangles by default).
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
plot.RutterGatsonisLCA <- function(x,
                                   symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                                   scale=0.02,
                                   size=c("eb","equal","sampsize","eb_revman","sampsize_revman"), 
                                   specrange=c(0.7,0.995),
                                   main="Diagnostic Test Accuracy Meta-Analysis", ...) {
  size    <- match.arg(size)
  symbol  <- match.arg(symbol)
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
  pointsXY(x=1-x$data$spec_eb, 
           y=x$data$sens_eb, 
           symbol = symbol, 
           scale = scale*0.5,
           cex.x = pct$sp,
           cex.y = pct$se,
           col="darkgray")
  symb <- switch(symbol, rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
  #points(x=XP$FPR,y=XP$sens,pch=0,col="darkgray",cex=2)
  # Add the ROC curve
  points(roc_points2, type="l", lwd=2,ann=FALSE)###
  # Add the legend 
  legend("bottomright", 
         bty ="n",
         legend = c(NA,
                    "HSROC curve",
                    "Empirical Bayes estimates"), 
         pch = c(NA,NA,symb), 
         lty = c(NA,1,NA), 
         lwd = c(NA,2,NA), 
         col = c(NA,"black","darkgray"))
  invisible(NULL)
}
