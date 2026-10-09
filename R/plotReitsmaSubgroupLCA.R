#' Plot Results from a Reitsma Subgroup LCA Model
#'
#' Produces a summary ROC plot for objects of class \code{"ReitsmaSubgroupLCA"}
#' obtained from \code{\link{fitReitsmaSubgroupLCA}}. Study-specific sensitivities and 
#' specificities are empirical Bayes estimates derived from the fitted latent 
#' class model. 
#'
#' @param x An object of class \code{"ReitsmaSubgroupLCA"}, as returned by
#'   \code{\link{fitReitsmaSubgroupLCA}}.
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
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{FALSE}.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param col Vector of colours used for subgroup-specific HSROC curves,
#'   study-level point estimates, summary points, confidence and prediction region. 
#'   If \code{NULL}, colours are generated automatically.
#' @param nudge_legend Numeric horizontal offset for the subgroup legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' 
#' @param conflevel Confidence level for the confidence region. Default is \code{0.95}.
#' @param predlevel Confidence level for the prediction region. Default is \code{0.95}.
#' @param connectstudies Whether the point estimates of two subgroups 
#'   within the same study should be connected. Defaults to \code{FALSE}.
#' @param studyCI Whether study-level confidence
#'   intervals should be displayed for the observed sensitivity and
#'   specificity estimates. Defaults to \code{FALSE}.
#'
#' @param studylabels Whether study identifiers should
#'   be displayed next to the study-level estimates. Defaults to
#'   \code{FALSE}.
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
#'   \item Study-level sensitivity and specificity estimates
#'   \item Summary (pooled) estimate
#'   \item confidence region around the summary point
#'   \item prediction region reflecting between-study variability
#'   \item Optional HSROC curve (if \code{HSROC = TRUE})
#' }
#'
#' Confidence and prediction regions are derived using the delta method
#' based on the estimated variance-covariance structure of the model.
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
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitReitsmaSubgroupLCA}}
#' @method plot ReitsmaSubgroupLCA
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.ReitsmaSubgroupLCA <- function(x, 
                                    symbol=NULL,
                                    scale=0.02, 
                                    size=c("eb","equal","sampsize","eb_revman","sampsize_revman"), 
                                    main="Diagnostic Test Accuracy Meta-Analysis",
                                    col=NULL,
                                    nudge_legend=-0.4,
                                    HSROC=FALSE,
                                    specrange=c(0.7,0.995),
                                    conflevel=0.95,
                                    predlevel=0.95,
                                    connectstudies=FALSE,
                                    studyCI=FALSE,
                                    studylabels=FALSE,
                                    ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  if (!is.numeric(predlevel) || length(predlevel) != 1L ||
      predlevel <= 0 || predlevel >= 1) {
    stop("predlevel must be a single number in (0, 1).")
  }
  if(connectstudies) {
    if(length(unique(x$data$subgroup)) != 2) {
      warning("'connectstudies=TRUE' is only recommended for two-subgroup comparisons." )
    }
  }
  size  <- match.arg(size)
  sub   <- levels(x$data$subgroup)
  nsub  <- length(sub)
  aes  <- validateSubgroupAesthetics(col=col,symbol=symbol,nsub=nsub)
  col  <- aes$col
  col2 <- aes$col2
  symbols2 <- aes$symbols2
  symb <- vector(mode="integer",length=nsub)
  # Calculations for percentage weights
  pct <- getWEIGHTSLCA(xdata=x$data,size=size)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  for (i in seq_along(sub)){
    sg    <- x$data$subgroup==sub[i]
    xspec <- x$data$spec_eb[sg]
    xsens <- x$data$sens_eb[sg]
    pctsp <- pct$sp[sg]
    pctse <- pct$se[sg]
    xstud <- x$data$study[sg]
    xdata <- x$data[sg,]
    ####
    symb[i] <- plotStudyLevelEstimates(xsens=xsens,xspec=xspec,xstud=xstud,xdata=xdata,
                                       symbol=symbols2[i],pctse=pctse,pctsp=pctsp,
                                       scale=scale,col=col2[i],
                                       conflevel=conflevel,studyCI=studyCI,
                                       studylabels=studylabels,LCA=TRUE)
  }
  # Add the ROC curve
  if(isTRUE(HSROC)){
    for(i in seq_along(sub)){
      roc_points2 <- getROCpoints(Lambda=x$RutterGatsonis_recovered[sub[i],"Lambda"],
                                  beta=x$RutterGatsonis_recovered[sub[i],"beta"],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=FALSE,col=col[i])
    }
  } ###
  # Add summary point
  for (i in seq_along(sub)){
    sg      <- sub[i]
    mu_A.sg <- paste0("mu_A.index.",sg)
    mu_B.sg <- paste0("mu_B.index.",sg)
    mean_point <- data.frame(1-x$sensspec[mu_B.sg,"Estimate"],
                             x$sensspec[mu_A.sg,"Estimate"])
    points(mean_point, col=col[i], cex=1.5, pch=15)
  }
  # Add confidence and prediction region
  for(i in seq_along(sub)){
    sg      <- sub[i]
    mu_A.sg <- paste0("mu_A.index.",sg)
    mu_B.sg <- paste0("mu_B.index.",sg) 
    s2_A.sg <- paste0("sigma2_A.index.",sg)
    s2_B.sg <- paste0("sigma2_B.index.",sg)
    s_AB.sg <- paste0("sigma_AB.index.",sg)
    muA     <- x$sdreport2[mu_A.sg,"Estimate"]
    muB     <- x$sdreport2[mu_B.sg,"Estimate"]
    seA     <- x$sdreport2[mu_A.sg,"Std. Error"]
    seB     <- x$sdreport2[mu_B.sg,"Std. Error"]
    covAB   <- x$vcov[mu_A.sg,mu_B.sg]
    varA    <- x$sdreport2[s2_A.sg,"Estimate"]
    varB    <- x$sdreport2[s2_B.sg,"Estimate"]
    sAB     <- x$sdreport2[s_AB.sg,"Estimate"]
    region  <- getConfPredRegion(muA=muA,muB=muB,
                                 seA=seA,seB=seB,covAB=covAB, # conf
                                 varA=varA,varB=varB,sAB=sAB, # pred
                                 nstudy=sum(x$data$subgroup == sg, na.rm = TRUE),
                                 conflevel=conflevel,
                                 predlevel=predlevel)
    lines(region$conf, lty=2, lwd=2, col=col2[i])
    lines(region$pred, lty=3, lwd=2, col=col2[i])
  }
  ## Connect studies
  if(isTRUE(connectstudies)){
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
  # Add the legend 
  conf_lab <- paste0(round(100 * conflevel), "% Confidence region")
  pred_lab <- paste0(round(100 * predlevel), "% Prediction region")
  if(isTRUE(HSROC)){
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "HSROC curve",
                      "Summary estimate",
                      conf_lab,
                      pred_lab,
                      "Empirical Bayes estimates"), 
           pch = c(NA,NA,15,NA,NA,symb[1]), 
           lty = c(NA,1,NA,2,3,NA), 
           lwd = c(NA,2,NA,2,2,NA), 
           col = c(NA,"black","black","black","black","darkgray"))}
  else{
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "Summary estimate",
                      conf_lab,
                      pred_lab,
                      "Empirical Bayes estimates"), 
           pch = c(NA,15,NA,NA,symb[1]), 
           lty = c(NA,NA,2,3,NA), 
           lwd = c(NA,NA,2,2,NA), 
           col = c(NA,"black","black","black","darkgray"))
  }
  legend("right",
         inset = c(nudge_legend, 0),
         legend = sub,
         col    = col,
         pch = symb,
         xpd = TRUE,
         cex = 1.2,
         bty = "n")
  invisible(NULL)
}
