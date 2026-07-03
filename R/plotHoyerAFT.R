#' Plot Results from a Hoyer Model
#'
#' Produces a hierarchical summary receiver operating characteristic
#' (HSROC) plot from a fitted Hoyer AFT model. The plot shows
#' study-level sensitivity and specificity estimates together with the
#' meta-analytic HSROC curve derived from the fitted model.
#'
#' @param x An object of class \code{"HoyerAFT"} as returned by
#'   \code{\link{fitHoyerAFT}}. Must contain:
#'   \describe{
#'     \item{original}{Original processed data including sensitivity
#'       (\code{sens}) and false positive rate (\code{fpr})}
#'     \item{rep2}{Summary of model parameters}
#'     \item{distcode}{Distribution code (1 = Weibull,
#'       2 = lognormal, 3 = loglogistic)}
#'   }
#'
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size. Default}
#'    \item{"sampsize"}{Size proportional to sample size}
#'    \item{"se"}{Size proportional to precision on the logit scale}
#'  }
#' @param thresholdrange A numeric vector of length 2 giving the range of
#'   threshold over which sensitivities and specificities are predicted
#'   If \code{NULL} (default), then the minimum and maximum thresholds
#'   from the data are used.
#'   
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#'   
#' @param ... Additional graphical arguments (currently unused).
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Study-specific sensitivity and false positive rate estimates
#'   \item Rectangular markers representing study observations
#'   \item Lines connecting thresholds within studies
#'   \item A meta-analytic HSROC curve based on the fitted AFT model
#' }
#'
#' The HSROC curve is constructed using the estimated model parameters
#' and depends on the specified distribution:
#' \itemize{
#'   \item Weibull
#'   \item Lognormal
#'   \item Loglogistic
#' }
#'
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#'
#' @seealso \code{\link{fitHoyerAFT}} \code{\link{fitHoyer}}
#' @importFrom graphics abline axis legend lines par points symbols title
#' @method plot HoyerAFT
#' @importFrom stats pnorm plogis
#' @export
plot.HoyerAFT <- function(x,scale=0.02, size=c("equal","sampsize","se"),
                          thresholdrange=NULL,
                          main="Diagnostic Test Accuracy Meta-Analysis",...) {
  size    <- match.arg(size)
  HH      <- x$data
  testdir <- unique(x$data$testdirection)
  if (length(testdir) != 1) stop("testdirection must be unique")
  
  op <- par(pty = "s")
  ### Plot coordinate system
  pct <- getWEIGHTS(HH,size)
  plot_SESPGRID(main=main)
  symbols(x=x$data$fpr,y=x$data$sens,rectangles=cbind(pct$sp,pct$se)*scale,inches=F,add=T,fg="darkgray")
  # Add lines
  studies <- unique(HH$study)
  for(i in seq_along(studies)) {
    # subset one study
    d <- HH[HH$study == studies[i], ]
    # order by threshold (or by FPR if you prefer)
    d <- d[order(d$threshold), ]
    # add connecting line
    lines(d$fpr, d$sens,
          col = "darkgray",
          lty = c(3,3),
          lwd = 1)
  }
  ### Plot meta-analytical summary ROC curve
  if(is.null(thresholdrange)){
    minth   <- min(HH$threshold)
    maxth   <- max(HH$threshold)
  } else {
    minth   <- thresholdrange[1]
    maxth   <- thresholdrange[2]
  }
  beta0   <- x$sdreport2["beta0","Estimate"]
  beta1   <- x$sdreport2["beta1","Estimate"]
  lambda0 <- x$sdreport2["lambda0","Estimate"]
  lambda1 <- x$sdreport2["lambda1","Estimate"]
  xx      <- exp(seq(log(minth),
                     log(maxth),
                     length.out = 1000))
  if(x$distcode==1 & testdir=="greater"){
    roc_points <- data.frame(fpr =exp(-(xx*exp(-(beta0)))**(1/lambda0)),
                             sens=exp(-(xx*exp(-(beta1)))**(1/lambda1)))}

  if(x$distcode==2 & testdir=="greater"){
    roc_points <- data.frame(fpr =1-pnorm((log(xx)-beta0)/lambda0),
                             sens=1-pnorm((log(xx)-beta1)/lambda1))}
  if(x$distcode==3 & testdir=="greater"){
    roc_points <- data.frame(fpr =plogis((beta0-log(xx))/lambda0),
                             sens=plogis((beta1-log(xx))/lambda1))}
  #############
  #############
  if(x$distcode==1 & testdir=="less"){
    roc_points <- data.frame(fpr =1-exp(-(xx*exp(-(beta0)))**(1/lambda0)),
                             sens=1-exp(-(xx*exp(-(beta1)))**(1/lambda1)))}
  if(x$distcode==2 & testdir=="less"){
    roc_points <- data.frame(fpr =pnorm((log(xx)-beta0)/lambda0),
                             sens=pnorm((log(xx)-beta1)/lambda1))}
  if(x$distcode==3 & testdir=="less"){
    roc_points <- data.frame(fpr =1-plogis((beta0-log(xx))/lambda0),
                             sens=1-plogis((beta1-log(xx))/lambda1))}
  ##########
  points(roc_points, type="l", lwd=2,ann=F)###
  # Add summary point
  # Add the legend
  legend("bottomright",
         bty ="n",
         legend = c(NA,"HSROC curve","Data"),
         pch = c(NA,NA,0),
         lty = c(NA,1,NA),
         lwd = c(NA,2,NA),
         col = c(NA,"black","darkgray"))
  on.exit(par(op))
}