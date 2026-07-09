#' Plot Results from a Reitsma Subgroup Model
#'
#' Produces a summary ROC plot for objects of class \code{"ReitsmaSubgroup"}
#' obtained from \code{\link{fitReitsmaSubgroup}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. Optionally, the
#' HSROC (hierarchical summary ROC) curve can be overlaid.
#'
#' @param x An object of class \code{"ReitsmaSubgroup"}, as returned by
#'   \code{\link{fitReitsmaSubgroup}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size. Default}
#'    \item{"sampsize"}{Size proportional to sample size}
#'    \item{"se"}{Size proportional to precision on the logit scale}
#'  }  
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{FALSE}.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param col Vector of colours used for subgroup-specific HSROC curves,
#'   study-level rectangles, summary points, confidence and prediction region. 
#'   If \code{NULL}, colours are generated automatically.
#' @param nudge_legend Numeric horizontal offset for the subgroup legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' 
#' @param conflevel Confidence level for the confidence region. Default is \code{0.95}.
#' @param predlevel Confidence level for the prediction region. Default is \code{0.95}.
#' @param connectstudies Whether the point estimates (rectangles) of two subgroups 
#'   within the same study should be connected. Defaults to \code{FALSE}.
#'
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown as rectangles, where the size reflects
#' approximate study weights derived from the Fisher information matrix.
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
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitReitsmaSubgroup}}
#' @method plot ReitsmaSubgroup
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.ReitsmaSubgroup <- function(x, scale=0.02, 
                                 size=c("equal","sampsize","se"), 
                                 main="Diagnostic Test Accuracy Meta-Analysis",
                                 col=NULL,
                                 nudge_legend=-0.4,
                                 HSROC=FALSE,
                                 specrange=c(0.7,0.995),
                                 conflevel=0.95,
                                 predlevel=0.95,
                                 connectstudies=FALSE,
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
  subs  <- levels(x$data$subgroup_safe)
  nsub  <- length(sub)
  if(!all(subs == make.names(sub))){
    stop("Object is corrupted. Please don't change object after running fitReitsmaSubgroup().")
  }
  if(!all(x$data$subgroup_safe == make.names(x$data$subgroup))){
    stop("Object is corrupted. Please don't change object after running fitReitsmaSubgroup().")
  }
  if(is.null(col)) col <- grDevices::rainbow(n=nsub)
  col2 <- grDevices::adjustcolor(col,alpha.f=0.6)
  # Calculations for percentage weights
  pct <- getWEIGHTS(xdata=x$data,size=size)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  for (i in seq_along(sub)){
    symbols(x=1-x$data$spec[x$data$subgroup==sub[i]],
            y=x$data$sens[x$data$subgroup==sub[i]],
            rectangles=cbind(pct$sp[x$data$subgroup==sub[i]],
                             pct$se[x$data$subgroup==sub[i]])*scale,
            inches=FALSE,
            add=TRUE,
            fg=col2[i])
  }
    # Add the ROC curve
  if(HSROC==TRUE){
    for(i in seq_along(sub)){
      roc_points2 <- getROCpoints(Lambda=x$RutterGatsonis_recovered[sub[i],"Lambda"],
                                  beta=x$RutterGatsonis_recovered[sub[i],"beta"],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=FALSE,col=col[i])
    }
  } ###
  # Add summary point
  for (i in seq_along(subs)){
    sg      <- subs[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    mean_point <- data.frame(1-x$sensspec[mu_B.sg,"Estimate"],
                             x$sensspec[mu_A.sg,"Estimate"])
    points(mean_point, col=col[i], cex=1.5, pch=15)
  }
  # Add confidence and prediction region
  for(i in seq_along(subs)){
    sg      <- subs[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    if(x$variances=="unequal"){
      s2_A.sg <- paste0("sigma2_A.",sg)
      s2_B.sg <- paste0("sigma2_B.",sg)
      s_AB.sg <- paste0("sigma_AB.",sg)
    }
    if(x$variances=="common"){
      s2_A.sg <- "sigma2_A.sens"
      s2_B.sg <- "sigma2_B.spec"
      s_AB.sg <- "sigma_AB"
    }
    muA     <- x$estimates_mu[mu_A.sg,"Estimate"]
    muB     <- x$estimates_mu[mu_B.sg,"Estimate"]
    seA     <- x$estimates_mu[mu_A.sg,"Std_Error"]
    seB     <- x$estimates_mu[mu_B.sg,"Std_Error"]
    covAB   <- x$vcov_mu[mu_A.sg,mu_B.sg]
    varA    <- x$estimates_mu[s2_A.sg,"Estimate"]
    varB    <- x$estimates_mu[s2_B.sg,"Estimate"]
    sAB     <- x$estimates_mu[s_AB.sg,"Estimate"]
    region  <- getConfPredRegion(muA=muA,muB=muB,
                                 seA=seA,seB=seB,covAB=covAB, # conf
                                 varA=varA,varB=varB,sAB=sAB, # pred
                                 nstudy=sum(x$data$subgroup_safe == sg, na.rm = TRUE),
                                 conflevel=conflevel,
                                 predlevel=predlevel)
    lines(region$conf, lty=2, lwd=2, col=col2[i])
    lines(region$pred, lty=3, lwd=2, col=col2[i])
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
  # Add the legend 
  conf_lab <- paste0(round(100 * conflevel), "% Confidence region")
  pred_lab <- paste0(round(100 * predlevel), "% Prediction region")
  if(HSROC==TRUE){
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "HSROC curve",
                      "Summary estimate",
                      conf_lab,
                      pred_lab,
                      "Data"), 
           pch = c(NA,NA,15,NA,NA,0), 
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
                      "Data"), 
           pch = c(NA,15,NA,NA,0), 
           lty = c(NA,NA,2,3,NA), 
           lwd = c(NA,NA,2,2,NA), 
           col = c(NA,"black","black","black","darkgray"))
  }
  legend("right",
         inset = c(nudge_legend, 0),
         legend = sub,
         col    = col,
         pch = 0,
         xpd = TRUE,
         cex = 1.2,
         bty = "n")
  invisible(NULL)
}
