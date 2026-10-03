#' Plot Results from a Rutter and Gatsonis Model
#'
#' Produces a summary ROC plot for objects of class \code{"RutterGatsonis"}
#' obtained from \code{\link{fitRutterGatsonis}}. The plot displays study-level
#' estimates of sensitivity and specificity and the
#' HSROC (hierarchical summary ROC).
#'
#' @param x An object of class \code{"RutterGatsonis"}, as returned by
#'   \code{\link{fitRutterGatsonis}}.
#' @param symbol Character string specifying the symbol used for study-level
#'   estimates. One of \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. The default is \code{"rectangle"}.
#' @param scale Optional numeric scaling factor controlling study-symbol
#'   dimensions. If \code{NULL} (default), \code{scale} is set to
#'   \code{0.02} for the original sizing methods and to \code{1} for
#'   RevMan-style sizing methods. For RevMan-style methods,
#'   \code{scale = 1} corresponds to the 100 percent point-scaling
#'   setting calibrated against a 540 by 540 SROC plotting panel.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size.}
#'    \item{"sampsize"}{
#'     Horizontal and vertical symbol dimensions are proportional to the
#'     relative numbers of non-diseased and diseased participants,
#'     respectively.}
#'    \item{"se"}{
#'       Horizontal and vertical symbol dimensions are proportional to
#'       normalized approximate inverse standard errors of observed logit
#'       specificity and logit sensitivity, respectively.}
#'    \item{"sampsize_revman"}{RevMan-style sample-size scaling. 
#'    Horizontal and vertical symbol dimensions are proportional to
#'    \eqn{4 + 0.64\sqrt{n_0}} and
#'    \eqn{4 + 0.64\sqrt{n_1}}, respectively, where \eqn{n_0} and
#'    \eqn{n_1} are the numbers of non-diseased and diseased participants.}
#'    \item{"se_revman"}{RevMan-style inverse-standard-error scaling. 
#'    Horizontal and vertical symbol dimensions are proportional to
#'    \eqn{4 + 1.5\sqrt{n_0\widehat{Sp}(1-\widehat{Sp})}} and
#'    \eqn{4 + 1.5\sqrt{n_1\widehat{Se}(1-\widehat{Se})}},
#'    respectively.}
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
#' Study-specific estimates are shown using the selected symbol (rectangles by default).
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates
#'   \item HSROC curve
#' }
#'
#' @references
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
#' @seealso \code{\link{fitRutterGatsonis}}
#' @method plot RutterGatsonis
#' @export
plot.RutterGatsonis <- function(x,
                                symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                                scale=NULL,
                                size=c("equal","sampsize","se","sampsize_revman","se_revman"), 
                                specrange=c(0.7,0.995),
                                main="Diagnostic Test Accuracy Meta-Analysis", ...) {
  size    <- match.arg(size)
  if (is.null(scale)) {
    scale <- if (size %in% c("sampsize_revman","se_revman")) {1} else {0.02} 
  }
  symbol  <- match.arg(symbol)
  warn_unestimable_sroc_points(x$data)
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
  pct <- getWEIGHTS(x$data,size)
  # Plot study level estimates 
  pointsXY(x=1-x$data$spec, 
           y=x$data$sens, 
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
                      "Data"), 
           pch = c(NA,NA,symb), 
           lty = c(NA,1,NA), 
           lwd = c(NA,2,NA), 
           col = c(NA,"black","darkgray"))
  invisible(NULL)
}
