#' Plot Results from a Rutter and Gatsonis Subgroup LCA Model
#'
#' Produces summary ROC plots for objects of class \code{"RutterGatsonisSubgroupLCA"}
#' obtained from \code{\link{fitRutterGatsonisSubgroupLCA}}. The plot displays
#' study-level estimates of sensitivity and specificity (empirical Bayes estimates), stratified by subgroup,
#' together with subgroup-specific HSROC (hierarchical summary ROC) curves.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroupLCA"}, as returned by
#'   \code{\link{fitRutterGatsonisSubgroupLCA}}.
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the subgroups. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally. Relative adjustments can be 
#' specified directly, for example, \code{scale = 0.02 * 0.3} produces marker 
#' dimensions equal to 30 percent of the default.
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
#'       RevMan-style graphical scaling applied to the inverse conditional
#'       posterior standard errors of the study-specific empirical Bayes
#'       logit-specificity and logit-sensitivity estimates. A common
#'       analysis-specific scaling constant is chosen so that the average
#'       larger dimension measure across studies equals eight before adding
#'       the fixed baseline dimension. This is an analogous graphical extension
#'       for latent-class models and is not an original Cochrane latent-class
#'       weighting method.}
#'     \item{\code{"sampsize_revman"}}{
#'       RevMan-style graphical scaling based on total study sample size.
#'       The same transformed total-sample-size measure is used horizontally
#'       and vertically.}
#'   }
#' @param col Vector of colours used for subgroup-specific HSROC curves
#'   and study-level point estimates. If \code{NULL}, colours are generated automatically.
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{TRUE}.
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
#' @param studyCI Whether study-level confidence
#'   intervals should be displayed for the observed sensitivity and
#'   specificity estimates. Defaults to \code{FALSE}.
#' @param conflevel Confidence level for the study-level confidence 
#'   intervals. Default is \code{0.95}.
#' @param studylabels Whether study identifiers should
#'   be displayed next to the study-level estimates. Defaults to
#'   \code{FALSE}.
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
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitRutterGatsonisSubgroupLCA}}
#' @method plot RutterGatsonisSubgroupLCA
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.RutterGatsonisSubgroupLCA <- function(x, 
                                           symbol=NULL,
                                           scale=0.02,
                                           size=c("eb","equal","sampsize","eb_revman","sampsize_revman"), 
                                           nudge_legend=-0.4,
                                           specrange=c(0.7,0.995),
                                           col=NULL,
                                           main="Diagnostic Test Accuracy Meta-Analysis",
                                           HSROC=TRUE,
                                           connectstudies=FALSE,
                                           studyCI=FALSE,
                                           conflevel=0.95,
                                           studylabels=FALSE,
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
  aes  <- validateSubgroupAesthetics(col=col,symbol=symbol,nsub=nsub)
  col  <- aes$col
  col2 <- aes$col2
  symbols2 <- aes$symbols2
  symb <- vector(mode="integer",length=nsub)
  ##
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")   # enlarge right margin
  plot_SESPGRID(main=main)
  # Data points
  pct <- getWEIGHTSLCA(xdata=x$data,size=size)
  for (i in seq_along(sub)){
    sg    <- x$data$subgroup==sub[i]
    xspec <- x$data$spec[sg]
    xsens <- x$data$sens[sg]
    pctsp <- pct$sp[sg]
    pctse <- pct$se[sg]
    xstud <- x$data$study[sg]
    ####
    pointsXY(x=1-xspec, 
             y=xsens, 
             symbol = symbols2[i], 
             scale = scale*0.5,
             cex.x = pctsp,
             cex.y = pctse,
             col=col2[i])
    symb[i] <- switch(symbols2[i], rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
    ###
    if(studyCI==TRUE) {
      forestci <- getForestSensSpecLCA(xdata=x$data[sg,],
                                       conflevel=conflevel)$XP[,c("sens",
                                                                  "spec",
                                                                  "Sens_LCI",
                                                                  "Sens_UCI",
                                                                  "Spec_LCI",
                                                                  "Spec_UCI")]
      drawStudyPointCI(x = 1 - forestci$spec,
                       y = forestci$sens,
                       sens_lower = forestci$Sens_LCI,
                       sens_upper = forestci$Sens_UCI,
                       spec_lower = forestci$Spec_LCI,
                       spec_upper = forestci$Spec_UCI,
                       symbol = symbols2[i],
                       col = col2[i],
                       scale = scale * 0.5,
                       cex.x = pctsp,
                       cex.y = pctse)}
    if(studylabels==TRUE){
      graphics::text(x = 1-xspec,
                     y = xsens,
                     cex = 1,
                     pos = 4,
                     col = col2[i],
                     labels = xstud)
    }
  }
  # Add the ROC curve
  if(HSROC==TRUE){
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
  }
  ## Connect studies
  if(connectstudies){
    for(st in unique(x$data$study)){
      tmp <- x$data[x$data$study == st, ]
      tmp <- tmp[order(tmp$subgroup), ]
      if(nrow(tmp) == 2) {
        lines(x = c(1-tmp$spec_eb[1],1-tmp$spec_eb[2]),
              y = c(tmp$sens_eb[1],tmp$sens_eb[2]),
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
  if(HSROC==TRUE){
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "HSROC curve",
                      "Empirical Bayes estimates"), 
           pch = c(NA,NA,symb[1]), 
           lty = c(NA,1,NA), 
           lwd = c(NA,2,NA), 
           col = c(NA,"black","darkgray"))
  }
  invisible(NULL)
}