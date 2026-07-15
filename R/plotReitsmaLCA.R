#' Plot Results from a Reitsma LCA Model
#'
#' Produces a summary ROC plot for objects of class \code{"ReitsmaLCA"}
#' obtained from \code{\link{fitReitsmaLCA}}. Study-specific sensitivities and 
#' specificities are empirical Bayes estimates derived from the fitted latent 
#' class model. 
#'
#' @param x An object of class \code{"ReitsmaLCA"}, as returned by
#'   \code{\link{fitReitsmaLCA}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"eb"}{Size proportional to the precision of the empirical Bayes estimates. Default.}
#'    \item{"equal"}{All studies shown with equal size}
#'    \item{"sampsize"}{Size proportional to sample size}
#'  }  
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{FALSE}.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#'
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' 
#' @param conflevel Confidence level for the confidence region. Default is \code{0.95}.
#' @param predlevel Confidence level for the prediction region. Default is \code{0.95}.
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
#'   \item Study-level sensitivity and specificity estimates (Empirical Bayes Estimates)
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
#' @seealso \code{\link{fitReitsmaLCA}}
#' @importFrom stats qlogis plogis predict qf
#' @method plot ReitsmaLCA
#' @export
plot.ReitsmaLCA <- function(x, scale=0.02, 
                            size=c("eb","equal","sampsize"), 
                            main="Diagnostic Test Accuracy Meta-Analysis",
                            HSROC=FALSE, 
                            specrange=c(0.7,0.995),
                            conflevel=0.95,
                            predlevel=0.95, ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  if (!is.numeric(predlevel) || length(predlevel) != 1L ||
      predlevel <= 0 || predlevel >= 1) {
    stop("predlevel must be a single number in (0, 1).")
  }
  size    <- match.arg(size)
  nstudy  <- nrow(x$data)
  # Confidence and prediction region
  muA     <- x$sdreport2["mu_A.index","Estimate"]
  muB     <- x$sdreport2["mu_B.index","Estimate"]
  seB     <- x$sdreport2["mu_B.index","Std. Error"]
  seA     <- x$sdreport2["mu_A.index","Std. Error"]
  covAB   <- x$vcov["mu_A.index","mu_B.index"]
  varA    <- x$sdreport2["sigma2_A.index","Estimate"]
  varB    <- x$sdreport2["sigma2_B.index","Estimate"]
  sAB     <- x$sdreport2["sigma_AB.index","Estimate"]
  region <- getConfPredRegion(muA=muA,muB=muB,
                              seA=seA,seB=seB,covAB=covAB, # conf
                              varA=varA,varB=varB,sAB=sAB, # pred
                              nstudy=nstudy,
                              conflevel=conflevel,
                              predlevel=predlevel)
  # Calculations for percentage weights
  pct <- getWEIGHTSLCA(x$data,size)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  symbols(x=1-x$data$spec_eb,y=x$data$sens_eb,rectangles=cbind(pct$sp,pct$se)*scale,inches=FALSE,add=TRUE,fg="darkgray")
  #points(x=XP$FPR,y=XP$sens,pch=0,col="darkgray",cex=2)
  # Add the ROC curve
  if(HSROC==TRUE){
    Lambda  <- x$RutterGatsonis_recovered$Lambda
    beta    <- x$RutterGatsonis_recovered$beta
    roc_points2 <- getROCpoints(Lambda,beta,specrange)
    points(roc_points2, type="l", lwd=2,ann=FALSE)
  } ###
  # Add summary point
  mean_point <- data.frame(1-x$sensspec["mu_B.index","Estimate"],
                           x$sensspec["mu_A.index","Estimate"])
  points(mean_point, col="black",cex=1.5, pch=15)
  # Add confidence and prediction region
  lines(region$conf, lty=2, lwd=2, col="black")
  lines(region$pred, lty=3, lwd=2, col="black")
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
                      "Empirical Bayes estimates"), 
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
                      "Empirical Bayes estimates"), 
           pch = c(NA,15,NA,NA,0), 
           lty = c(NA,NA,2,3,NA), 
           lwd = c(NA,NA,2,2,NA), 
           col = c(NA,"black","black","black","darkgray"))
  }
  invisible(NULL)
}
