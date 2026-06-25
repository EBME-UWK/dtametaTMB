#' Plot Results from a Rutter and Gatsonis Diagnostic Test Accuracy Model
#'
#' Produces a summary ROC plot for objects of class \code{"RutterGatsonis"}
#' obtained from \code{\link{fitRutterGatsonis}}. The plot displays study-level
#' estimates of sensitivity and specificity and the
#' HSROC (hierarchical summary ROC).
#'
#' @param x An object of class \code{"RutterGatsonis"}, as returned by
#'   \code{\link{fitRutterGatsonis}}.
#' @param scale A numeric scaling factor controlling the size of the
#'   rectangles representing study weights. Default is \code{0.02}.
#' @param size Character string controlling study weight display:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size. Default}
#'    \item{"sampsize"}{Size proportional to sample size}
#'    \item{"se"}{Size proportional to precision on the logit scale}
#'  }
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 − false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown as rectangles.
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates
#'   \item HSROC curve
#' }
#'
#' @return
#' This function is called for its side effect of producing a plot and
#' returns \code{invisible(NULL)}.
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
#' @seealso \code{\link{fitRutterGatsonis}}
#' @importFrom stats qlogis plogis 
#' @method plot RutterGatsonis
#' @export

plot.RutterGatsonis <- function(x, scale=0.02,size=c("equal","sampsize","se"), ...) {
  size    <- match.arg(size)
  Lambda  <- x$sdreport2["Lambda", "Estimate"]
  beta    <- x$sdreport2["beta","Estimate"]
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
  title(main="Diagnostic Accuracy Meta-Analysis", xlab="Specificity", ylab="Sensitivity")
  # Plot study level estimates 
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
  symbols(x=1-x$data$spec,y=x$data$sens,rectangles=cbind(pctsp,pctse)*scale,inches=F,add=T,fg="darkgray")
  #points(x=XP$FPR,y=XP$sens,pch=0,col="darkgray",cex=2)
  # Add the ROC curve
  points(roc_points2, type="l", lwd=2,ann=F)###
  # Add the legend 
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
