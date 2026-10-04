#' Plot Results from a Reitsma Model
#'
#' Produces a summary ROC plot for objects of class \code{"Reitsma"}
#' obtained from \code{\link{fitReitsma}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. Optionally, the
#' HSROC (hierarchical summary ROC) curve can be overlaid.
#'
#' @param x An object of class \code{"Reitsma"}, as returned by
#'   \code{\link{fitReitsma}}.
#' @param symbol Character string specifying the symbol used for study-level
#'   estimates. One of \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. The default is \code{"rectangle"}.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"fisher"}{Size proportional to a decomposition of Fisher's Information matrix. Default.}
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
#' Study-specific estimates are shown using the selected symbol (rectangles by default), where the size reflects
#' approximate study weights derived from the Fisher information matrix (default).
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
#' Riley, R. D., Ensor, J., Jackson, D., & Burke, D. L. (2018).
#' Deriving percentage study weights in multi-parameter meta-analysis models:
#' with application to meta-regression, network meta-analysis and one-stage
#' individual participant data models.
#' \emph{Statistical Methods in Medical Research}, 27(10), 2885--2905.
#' \doi{10.1177/0962280216688033}
#'
#' @return
#' No return value. Called for its side effect of producing a plot.
#' @seealso \code{\link{fitReitsma}}
#' @importFrom stats qlogis plogis predict qf
#' @method plot Reitsma
#' @export
plot.Reitsma <- function(x, 
                         symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                         scale=0.02, 
                         size=c("fisher","equal","sampsize","se","sampsize_revman","se_revman"), 
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
  symbol  <- match.arg(symbol)
  size    <- match.arg(size)
  nstudy  <- nrow(x$data)
  warn_unestimable_sroc_points(x$data)
  # Confidence and prediction region
  muA     <- x$estimates["mu_A.sens","Estimate"]
  muB     <- x$estimates["mu_B.spec","Estimate"]
  seB     <- x$estimates["mu_B.spec","Std_Error"]
  seA     <- x$estimates["mu_A.sens","Std_Error"]
  covAB   <- x$vcov["mu_A.sens","mu_B.spec"]
  varA    <- x$estimates["sigma2_A.sens","Estimate"]
  varB    <- x$estimates["sigma2_B.spec","Estimate"]
  sAB     <- x$estimates["sigma_AB","Estimate"]
  region <- getConfPredRegion(muA=muA,muB=muB,
                              seA=seA,seB=seB,covAB=covAB, # conf
                              varA=varA,varB=varB,sAB=sAB, # pred
                              nstudy=nstudy,
                              conflevel=conflevel,
                              predlevel=predlevel)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  ###
  # Calculations for percentage weights
  if(size=="fisher"){
    X <- x$data
    X$n1    <- X$TP+X$FN
    X$n0    <- X$FP+X$TN
    X$true1 <- X$TP
    X$true0 <- X$TN 
    X$recordid <- seq_len(nrow(X))
    Y_pw <- reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1", "true0")), 
                    timevar="sens", times=c(1,0), v.names=c("n","true")) 
    ##
    Y_pw = Y_pw[order(Y_pw$id),]
    Y_pw$spec <- 1-Y_pw$sens
    X_pw <- cbind(Y_pw$sens, Y_pw$spec)
    XT_pw <- t(X_pw)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    invn[!is.finite(invn)] <- 100000000000 
    A <- diag(invn)
    p_pw <- stats::predict(x$glmmTMB, type="response")
    var_pw <- p_pw*(1-p_pw)
    B <- diag(var_pw)
    #G_one <- matrix(c(varA,varAB,varAB,varB),2,2)
    G_one <- matrix(c(varA,sAB,sAB,varB),2,2)
    G <- kronecker(diag(nstudy), G_one)
    #inverse of B (required later on)
    BI <- solve(B)
    # Create variance matrix for observations
    V <- (Z %*% G %*% t(Z)) + (A %*% BI)
    # invert the variance matrix
    invV <- solve(V)
    # derive the fishers information matrix
    fish <- XT_pw %*% invV %*% X_pw
    # invert Fishers information to obtain Var Beta hat
    varb   <- solve(fish)
    pct <- data.frame(sp = numeric(nstudy),
                      se = numeric(nstudy))
    # Get weights  
    for (i in seq_len(nstudy)){
      DM <- V
      DM[(i*2)-1, (i*2)-1] <- 100000000000
      DM[(i*2)-1, (i*2)] <- 0
      DM[(i*2), (i*2)-1] <- 0
      DM[(i*2), (i*2)] <- 100000000000
      invDM <- solve(DM)
      fishD <- XT_pw %*% invDM %*% X_pw
      fishI <- fish - fishD
      weight <- varb %*% fishI %*% varb
      pct$se[i] <- 100*(weight[1,1]/varb[1,1])
      pct$sp[i] <- 100*(weight[2,2]/varb[2,2])
    }
    pct$se <- pmax(0,pct$se)
    pct$sp <- pmax(0,pct$sp)
  } else {
  pct <- getWEIGHTS(xdata=x$data,size=size)
  }
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
  if(HSROC==TRUE){
    Lambda  <- x$RutterGatsonis_recovered$Lambda
    beta    <- x$RutterGatsonis_recovered$beta
    roc_points2 <- getROCpoints(Lambda,beta,specrange)
    points(roc_points2, type="l", lwd=2,ann=FALSE)
    } ###
  # Add summary point
  mean_point <- data.frame(1-x$sensspec["spec","Estimate"],
                           x$sensspec["sens","Estimate"])
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
                      "Data"), 
          pch = c(NA,NA,15,NA,NA,symb), 
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
           pch = c(NA,15,NA,NA,symb), 
           lty = c(NA,NA,2,3,NA), 
           lwd = c(NA,NA,2,2,NA), 
           col = c(NA,"black","black","black","darkgray"))
  }
  invisible(NULL)
}
