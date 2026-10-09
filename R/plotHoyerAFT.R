#' Plot Results from a HoyerAFT Model
#'
#' Produces a hierarchical summary receiver operating characteristic
#' (HSROC) plot from a fitted HoyerAFT model. The plot shows
#' study-level sensitivity and specificity estimates together with the
#' meta-analytic HSROC curve derived from the fitted model.
#'
#' @param x An object of class \code{"HoyerAFT"} as returned by
#'   \code{\link{fitHoyerAFT}}. Must contain:
#'   \describe{
#'     \item{data}{Original processed data including sensitivity
#'       (\code{sens}) and false positive rate (\code{fpr})}
#'     \item{sdreport2}{Summary of model parameters}
#'     \item{distcode}{Distribution code (1 = Weibull,
#'       2 = lognormal, 3 = loglogistic)}
#'   }
#' @param symbol Character string specifying the symbol used for study-level
#'   estimates. One of \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. The default is \code{"rectangle"}.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally. Relative adjustments can be 
#' specified directly, for example, \code{scale = 0.02 * 0.3} produces marker 
#' dimensions equal to 30 percent of the default.
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
#'    \item{"sampsize_revman"}{
#'       RevMan-style sample-size scaling. Horizontal and vertical symbol
#'       dimensions are based on the square roots of the numbers of
#'       non-diseased and diseased participants, respectively. A common
#'       analysis-specific scaling constant is chosen so that the average
#'       larger dimension measure across studies equals eight before adding
#'       the fixed baseline dimension.}
#'    \item{"se_revman"}{
#'       RevMan-style inverse-standard-error scaling. Horizontal and
#'       vertical symbol dimensions are based on approximate inverse standard
#'       errors of observed logit specificity and logit sensitivity,
#'       respectively. A common analysis-specific scaling constant is chosen
#'       so that the average larger dimension measure across studies equals
#'       eight before adding the fixed baseline dimension.}
#'  } 
#'  
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{TRUE}.
#' @param thresholdrange A numeric vector of length 2 giving the range of
#'   threshold over which sensitivities and specificities are predicted.
#'   If \code{NULL} (default), then the minimum and maximum thresholds
#'   from the data are used.
#'   
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#'   
#' @param studyCI Whether study-level confidence
#'   intervals should be displayed for the observed sensitivity and
#'   specificity estimates. Defaults to \code{FALSE}.
#' @param conflevel Confidence level for the study-level confidence 
#'   intervals. Default is \code{0.95}.
#'
#' @param studylabels Whether study identifiers should
#'   be displayed next to the study-level estimates. Defaults to
#'   \code{FALSE}.
#'   
#' @param ... Additional graphical arguments (currently unused).
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Study-specific sensitivity and false positive rate estimates
#'   \item Symbols representing study observations, with dimensions determined by the selected graphical sizing method.
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
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitHoyerAFT}} \code{\link{fitHoyer}}
#' @importFrom graphics abline axis legend lines par points title
#' @importFrom stats plogis pnorm
#' @method plot HoyerAFT
#' @export
plot.HoyerAFT <- function(x,
                          symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                          scale=0.02, 
                          size=c("equal","sampsize","se","sampsize_revman","se_revman"),
                          HSROC=TRUE,
                          thresholdrange=NULL,
                          main="Diagnostic Test Accuracy Meta-Analysis",
                          studyCI=FALSE,
                          conflevel=0.95,
                          studylabels=FALSE,
                          ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  size    <- match.arg(size)
  symbol  <- match.arg(symbol)
  HH      <- x$data
  testdir <- unique(x$data$testdirection)
  if (length(testdir) != 1) stop("testdirection must be unique")
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(pty = "s")
  ### Plot coordinate system
  pct <- getWEIGHTS(HH,size)
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  # Plot study level estimates 
  xspec <- x$data$spec
  xsens <- x$data$sens
  pctsp <- pct$sp
  pctse <- pct$se
  xstud <- x$data$study
  xdata <- x$data
  symb <- plotStudyLevelEstimates(xsens=xsens,xspec=xspec,xstud=xstud,xdata=xdata,
                                  symbol=symbol,pctse=pctse,pctsp=pctsp,
                                  scale=scale,col="darkgray",
                                  conflevel=conflevel,studyCI=studyCI,
                                  studylabels=studylabels,LCA=FALSE)
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
    roc_points <- data.frame(fpr =1-stats::pnorm((log(xx)-beta0)/lambda0),
                             sens=1-stats::pnorm((log(xx)-beta1)/lambda1))}
  if(x$distcode==3 & testdir=="greater"){
    roc_points <- data.frame(fpr =stats::plogis((beta0-log(xx))/lambda0),
                             sens=stats::plogis((beta1-log(xx))/lambda1))}
  #############
  #############
  if(x$distcode==1 & testdir=="less"){
    roc_points <- data.frame(fpr =1-exp(-(xx*exp(-(beta0)))**(1/lambda0)),
                             sens=1-exp(-(xx*exp(-(beta1)))**(1/lambda1)))}
  if(x$distcode==2 & testdir=="less"){
    roc_points <- data.frame(fpr =stats::pnorm((log(xx)-beta0)/lambda0),
                             sens=stats::pnorm((log(xx)-beta1)/lambda1))}
  if(x$distcode==3 & testdir=="less"){
    roc_points <- data.frame(fpr =1-stats::plogis((beta0-log(xx))/lambda0),
                             sens=1-stats::plogis((beta1-log(xx))/lambda1))}
  ##########
  if(isTRUE(HSROC)){
    points(roc_points, type="l", lwd=2,ann=FALSE)###
  # Add the legend
  legend("bottomright",
         bty ="n",
         legend = c(NA,"HSROC curve","Data"),
         pch = c(NA,NA,symb),
         lty = c(NA,1,NA),
         lwd = c(NA,2,NA),
         col = c(NA,"black","darkgray"))
  }
  invisible(NULL)
}