#' Plot Results from a Rutter and Gatsonis Subgroup Model
#'
#' Produces summary ROC plots for objects of class \code{"RutterGatsonisSubgroup"}
#' obtained from \code{\link{fitRutterGatsonisSubgroup}}. The plot displays
#' study-level estimates of sensitivity and specificity, stratified by subgroup,
#' together with subgroup-specific HSROC (hierarchical summary ROC) curves.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroup"}, as returned by
#'   \code{\link{fitRutterGatsonisSubgroup}}.
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the subgroups. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally.
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
#' @param col Vector of colours used for subgroup-specific HSROC curves
#'   and study-level point estimates. If \code{NULL}, colours are generated automatically.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param nudge_legend Numeric horizontal offset for the subgroup legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param connectstudies Whether the point estimates of two subgroups 
#'   within the same study should be connected. Defaults to \code{FALSE}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#'   
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown with subgroup-specific symbols and colours.
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
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitRutterGatsonisSubgroup}}
#' @method plot RutterGatsonisSubgroup
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.RutterGatsonisSubgroup <- function(x, 
                                        symbol=NULL,
                                        scale=0.02,
                                        size=c("equal","sampsize","se","sampsize_revman","se_revman"), 
                                        nudge_legend=-0.4,
                                        specrange=c(0.7,0.995),
                                        col=NULL,
                                        main="Diagnostic Test Accuracy Meta-Analysis",
                                        connectstudies=FALSE,
                                        ...){
   if(connectstudies) {
     if(length(unique(x$data$subgroup)) != 2) {
       warning("'connectstudies=TRUE' is only recommended for two-subgroup comparisons." )
     }
   }
   size <- match.arg(size)
   sub  <- x$subgroups
   nsub <- length(sub)
   nstudy <- nrow(x$data)
   warn_unestimable_sroc_points(x$data)
   if(is.null(col)) col <- grDevices::rainbow(n=nsub)
   col2 <- grDevices::adjustcolor(col,alpha.f=0.6)
   # symbols
   if(is.null(symbol)) {
     available_symbols <- c("rectangle","ellipse","diamond","triangle","cross","plus","star")
     symbols2 <- rep(available_symbols,length.out=nsub)
   } else {
     symbols2 <- symbol
   }
   symb <- vector(mode="integer",length=nsub)
   ##
   oldpar <- par(no.readonly = TRUE)
   on.exit(par(oldpar))
   par(mar = c(5, 4, 4, 10),
       pty="s")   # enlarge right margin
   plot_SESPGRID(main=main)
   # Data points
   pct <- getWEIGHTS(xdata=x$data,size=size)
   for (i in seq_along(sub)){
     pointsXY(x=1-x$data$spec[x$data$subgroup==sub[i]], 
              y=x$data$sens[x$data$subgroup==sub[i]], 
              symbol = symbols2[i], 
              scale = scale*0.5,
              cex.x = pct$sp[x$data$subgroup==sub[i]],
              cex.y = pct$se[x$data$subgroup==sub[i]],
              col=col2[i])
     symb[i] <- switch(symbols2[i], rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
   }
   ###
   lamb <- paste0("Lambda_",sub)
   bet  <- paste0("beta_",sub)
  
   Lambda <- x$sdreport2[lamb,]
   beta   <- x$sdreport2[bet,]
   for(i in seq_len(nsub)){
      roc_points2 <- getROCpoints(Lambda[i,"Estimate"],
                                  beta[i,"Estimate"],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=FALSE,col=col[i])
   }
   ## Connect studies
   if(connectstudies){
     for(st in unique(x$data$study)){
       tmp <- x$data[x$data$study == st, ]
       tmp <- tmp[order(tmp$subgroup), ]
       if(nrow(tmp) == 2) {
         lines(x = c(1-tmp$spec[1],1-tmp$spec[2]),
               y = c(tmp$sens[1],tmp$sens[2]),
               type="l",col = "grey70",lwd =1)
       }
     }
   }
   ###
   legend("right",
          inset = c(nudge_legend, 0),
          legend = sub,
          col    = col,
          pch = symb,
          xpd = TRUE,
          cex = 1.2,
          bty = "n")

   legend("bottomright", 
          bty ="n",
          legend = c(NA,
                     "HSROC curve",
                     "Data"), 
          pch = c(NA,NA,symb[1]), 
          lty = c(NA,1,NA), 
          lwd = c(NA,2,NA), 
          col = c(NA,"black","darkgray"))
   invisible(NULL)
}