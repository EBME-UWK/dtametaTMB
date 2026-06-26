#' Plot Results from a Reitsma Diagnostic Test Accuracy Model
#'
#' Produces a summary ROC plot for objects of class \code{"Reitsma"}
#' obtained from \code{\link{fitReitsma}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. Optionally, the
#' HSROC (hierarchical summary ROC) curve can be overlaid.
#'
#' @param x An object of class \code{"Reitsma"}, as returned by
#'   \code{\link{fitReitsma}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"fisher"}{Size proportional to a decomposition of Fisher's Information matrix. Default}
#'    \item{"equal"}{All studies shown with equal size}
#'    \item{"sampsize"}{Size proportional to sample size}
#'    \item{"se"}{Size proportional to precision on the logit scale}
#'  }  
#' @param HSROC Logical; if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{FALSE}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' 
#' @param conflevel Numeric. Confidence level for the confidence region
#'   of the summary estimate. Must be between 0 and 1. Default is \code{0.95}.
#' @param predlevel Numeric. Confidence level for the prediction region.
#'   Must be between 0 and 1. Default is \code{0.95}.
#'
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 − false positive rate on a reversed axis).
#'
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
#' Freeman, S. C., Kerby, C. R., Patel, A., Cooper, N. J.,
#' Quinn, T., & Sutton, A. J. (2019).
#' Development of an interactive web-based tool to conduct
#' and interrogate meta-analysis of diagnostic test accuracy studies:
#' MetaDTA.
#' \emph{BMC Medical Research Methodology}, 19, 81.
#' \doi{10.1186/s12874-019-0724-x}
#'
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
#' This function is called for its side effect of producing a plot and
#' returns \code{invisible(NULL)}.
#'
#' @seealso \code{\link{fitReitsma}}
#' @importFrom stats qlogis plogis predict qf
#' @method plot Reitsma
#' @export


plot.Reitsma <- function(x, scale=0.02, 
                            size=c("fisher","equal","sampsize","se"), 
                            main="Diagnostic Test Accuracy Meta-Analysis",
                            HSROC=FALSE, 
                            conflevel=0.95,
                            predlevel=0.95, ...) {
  size    <- match.arg(size)
  nstudy  <- nrow(x$data)
  Lambda  <- x$RutterGatsonis_recovered$Lambda
  beta    <- x$RutterGatsonis_recovered$beta
  roc_points <- c()
  for (i in seq(from=0.001, to=0.999, by=0.01)){
    Sp_i  <- i
    Fpr_i <- 1-Sp_i
    LSp_i <- qlogis(Sp_i)
    LSe_i <- Lambda*exp(-beta/2) - exp(-beta)*LSp_i
    Se_i  <- plogis(LSe_i)
    roc_i <- data.frame(FPR=Fpr_i, Sen=Se_i)
    roc_points<-rbind(roc_points,roc_i)
  }
  minSens <- min(x$data$sens)
  maxSens <- max(x$data$sens)
  minFPR  <- min(1-x$data$spec)
  maxFPR  <- max(1-x$data$spec)
  # Create new data frame which restricts roc_points to being between min and max values
  roc_points2 <- 
     roc_points[
        roc_points$FPR < maxFPR &
        roc_points$FPR > minFPR &
        roc_points$Sen < maxSens &
        roc_points$Sen > minSens, ]
  
  # Confidence and prediction region
  muA     <- x$estimates["mu_A.sens",]$Estimate
  muB     <- x$estimates["mu_B.spec",]$Estimate
  seB     <- x$estimates["mu_B.spec",]$Std_Error
  seA     <- x$estimates["mu_A.sens",]$Std_Error
  covAB   <- x$vcov[1,2]
  r       <- covAB / (seA*seB)
  varA    <- x$estimates["sigma2_A.sens",]$Estimate
  varB    <- x$estimates["sigma2_B.spec",]$Estimate
  sAB     <- x$estimates["sigma_AB",]$Estimate
  sepredA <- sqrt(varA + seA**2)
  sepredB <- sqrt(varB + seB**2)
  rpredAB <- (sAB + covAB) / (sepredA*sepredB)
  f_conf  <- qf(conflevel, df1 = 2, df2 = nstudy - 2)
  f_pred  <- qf(predlevel, df1 = 2, df2 = nstudy - 2)
  croot_conf <- sqrt(2 * f_conf)
  croot_pred <- sqrt(2 * f_pred)
  
  conf_region <- c()
  pred_region <- c()
  # Confidence region
  for (i in seq(0, 2*pi, length.out=361)){
    confA  <- muA + (seA*croot_conf*cos(i))
    confB  <- muB + (seB*croot_conf*cos(i + acos(r)))
    confsens <- plogis(confA)
    confspec <- plogis(confB)
    conf_i <- data.frame(X=1-confspec, Y=confsens)
    #conf_i <- logit(conf_i)
    conf_region<-rbind(conf_region, conf_i)
  }
  for (i in seq(0, 2*pi, length.out=361)){
    predA <- muA + (sepredA*croot_pred*cos(i))
    predB <- muB + (sepredB*croot_pred*cos(i + acos(rpredAB)))
    predsens <- plogis(predA)
    predspec <- plogis(predB)
    pred_i <- data.frame(X=1-predspec, Y=predsens)
    #pred_i <- logit(pred_i)
    pred_region<-rbind(pred_region, pred_i)
  }
  # Calculations for percentage weights
  if(size=="fisher"){
    X <- x$data
    X$n1    <- X$TP+X$FN
    X$n0    <- X$FP+X$TN
    X$true1 <- X$TP
    X$true0 <- X$TN 
    X$recordid <- 1:nrow(X)
    Y_pw <- reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1", "true0")), 
                    timevar="sens", times=c(1,0), v.names=c("n","true")) 
    ##
    Y_pw = Y_pw[order(Y_pw$id),]
    Y_pw$spec <- 1-Y_pw$sens
    X_pw <- cbind(Y_pw$sens, Y_pw$spec)
    XT_pw <- t(X_pw)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    A <- diag(invn)
    p_pw <- stats::predict(x$glmmTMB, type="response")
    var_pw <- p_pw*(1-p_pw)
    B <- diag(var_pw)
    G_one <- matrix(c(varA,covAB,covAB,varB),2,2)
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
    varb  <- solve(fish)
    pctse <- vector(mode="numeric", length =nstudy)
    pctsp <- vector(mode="numeric", length =nstudy)
    # Get weights  
    for (i in 1:nstudy){
      DM <- V
      # DM2 <- diag(rep(100000000000),72)
      # DM1 <- V[c(1:2),c(1:2)]
      # VD1 <- Matrix::bdiag(DM1,DM2)
      # invDM <- solve(VD1)
      DM[(i*2)-1, (i*2)-1] <- 100000000000
      DM[(i*2)-1, (i*2)] <- 0
      DM[(i*2), (i*2)-1] <- 0
      DM[(i*2), (i*2)] <- 100000000000
      invDM <- solve(DM)
      fishD <- XT_pw %*% invDM %*% X_pw
      fishI <- fish - fishD
      weight <- varb %*% fishI %*% varb
      pctse[i] <- 100*(weight[1,1]/varb[1,1])
      pctsp[i] <- 100*(weight[2,2]/varb[2,2])
    }
  }
  if(size=="equal"){
    pctse <- rep(1,nrow(x$data))
    pctsp <- rep(1,nrow(x$data))
  }
  if(size=="sampsize"){
    pctse <- x$data$n1 / sum(x$data$n1)*100
    pctsp <- x$data$n0 / sum(x$data$n0)*100
  }
  if(size=="se"){
    sem1  <- x$data$sens*(1-x$data$sens)*x$data$n1 # inverse logit variance
    spm1  <- x$data$spec*(1-x$data$spec)*x$data$n0 # inverse logit variance
    pctse <- sqrt(sem1) / sum(sqrt(sem1))*100
    pctsp <- sqrt(spm1) / sum(sqrt(spm1))*100
  }
  ####
  op <- par(pty = "s")
  plot(1,1, ylim=c(0,1), xlim=c(0,1), xaxt = "n", yaxt="n",
       ann=F, pch=20, col="white",las=1,asp=1)
  axis( side = 1,                          # 1 = bottom axis
        at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  # positions of ticks
        labels = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0))  # custom labels
  axis( side = 2,
        at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),las=1)
  par(new=TRUE) 
  abline(v=(seq(0,1,0.2)), col="lightgray", lty="dotted")
  abline(h=(seq(0,1,0.2)), col="lightgray", lty="dotted")
  lines(c(0,1),c(0,1),col="lightgray",lty="dotted")
  # Add titles
  title(main=main, xlab="Specificity", ylab="Sensitivity")
  # Plot study level estimates 
  symbols(x=1-x$data$spec,y=x$data$sens,rectangles=cbind(pctsp,pctse)*scale,inches=F,add=T,fg="darkgray")
  #points(x=XP$FPR,y=XP$sens,pch=0,col="darkgray",cex=2)
  # Add the ROC curve
  if(HSROC==TRUE){points(roc_points2, type="l", lwd=2,ann=F)} ###
  # Add summary point
  mean_point <- data.frame(1-x$sensspec["spec",]$Estimate,
                           x$sensspec["sens",]$Estimate)
  points(mean_point, col="black",cex=1.5, pch=15)
  # Add confidence and prediction region
  lines(conf_region, lty=2, lwd=2, col="black")
  lines(pred_region, lty=3, lwd=2, col="black")
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
  on.exit(par(op))
}
