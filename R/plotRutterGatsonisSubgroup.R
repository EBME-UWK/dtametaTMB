#' Plot Results from a Rutter and Gatsonis Subgroup Model
#'
#' Produces summary ROC plots for objects of class \code{"RutterGatsonisSubgroup"}
#' obtained from \code{\link{fitRutterGatsonisSubgroup}}. The plot displays
#' study-level estimates of sensitivity and specificity, stratified by subgroup,
#' together with subgroup-specific HSROC (hierarchical summary ROC) curves.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroup"}, as returned by
#'   \code{\link{fitRutterGatsonisSubgroup}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'   \describe{
#'     \item{"equal"}{All studies shown with equal size. Default.}
#'     \item{"sampsize"}{Size proportional to sample size.}
#'     \item{"se"}{Size proportional to precision on the logit scale.}
#'   }
#' @param col Vector of colours used for subgroup-specific HSROC curves
#'   and study-level rectangles. If \code{NULL}, colours are generated automatically.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param nudge_legend Numeric value controlling the horizontal position of the
#'   subgroup legend relative to the right side of the plotting region.
#'   More negative values move the legend further to the right (outside the plot),
#'   whereas values closer to zero move it closer to the plotting area.
#'   Default is \code{-0.4}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown as rectangles, with subgroup-specific colours.
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates by subgroup
#'   \item Subgroup-specific HSROC curves
#'   \item A legend identifying the subgroups
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
#' @seealso \code{\link{fitRutterGatsonisSubgroup}}
#' @method plot RutterGatsonisSubgroup
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.RutterGatsonisSubgroup <- function(x, 
                                        scale=0.02,
                                        size=c("equal","sampsize","se"), 
                                        nudge_legend=-0.4,
                                        specrange=c(0.7,0.995),
                                        col=NULL,
                                        main="Diagnostic Test Accuracy Meta-Analysis",
                                        ...){
   size <- match.arg(size)
   sub  <- x$subgroups
   nsub <- length(sub)
   nstudy <- nrow(x$data)
   if(is.null(col)) col <- rainbow(n=nsub)
   col2 <- adjustcolor(col,alpha.f=0.6)
   
   op <- par(mar = c(5, 4, 4, 10),
             pty="s")   # enlarge right margin
   plot_SESPGRID(main=main)
   # Data points
   pct <- getWEIGHTS(xdata=x$data,size=size)
   for (i in seq_len(nsub)){
     symbols(x=1-x$data$spec[x$data$subgroup==sub[i]],
             y=x$data$sens[x$data$subgroup==sub[i]],
             rectangles=cbind(pct$sp[x$data$subgroup==sub[i]],
                              pct$se[x$data$subgroup==sub[i]])*scale,
             inches=F,
             add=T,
             fg=col2[i])
   }
   
   lamb <- paste0("Lambda_",sub)
   bet  <- paste0("beta_",sub)
  
   Lambda <- x$sdreport2[lamb,]
   beta   <- x$sdreport2[bet,]
   for(i in seq_len(nsub)){
      roc_points2 <- getROCpoints(Lambda[i],
                                  beta[i],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=F,col=col[i])
   }
   
   legend("right",
          inset = c(nudge_legend, 0),
          legend = sub,
          col    = col,
          pch = 0,
          xpd = TRUE,
          cex = 1.2,
          bty = "n")

   legend("bottomright", 
          bty ="n",
          legend = c(NA,
                     "HSROC curve",
                     "Data"), 
          pch = c(NA,NA,0), 
          lty = c(NA,1,NA), 
          lwd = c(NA,2,NA), 
          col = c(NA,"black","darkgray"))
   
    on.exit(par(op))
}