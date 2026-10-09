#' Plot Results from a Rutter and Gatsonis Model
#'
#' Produces a summary ROC plot for objects of class \code{"RutterGatsonis"}
#' obtained from \code{\link{fitRutterGatsonis}}. The plot displays study-level
#' estimates of sensitivity and specificity and the
#' HSROC (hierarchical summary ROC).
#'
#' @param x An object of class \code{"RutterGatsonis"}, as returned by
#'   \code{\link{fitRutterGatsonis}}.
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
#'    \item{"fisher"}{Size proportional to a decomposition of Fisher's Information matrix 
#'    with respect to the recovered pooled logit sensitivity and specificity. Default.}
#'    \item{"fisher_revman"}{
#'       Horizontal and vertical symbol dimensions are based on the same
#'       per-study percentage contributions to the logit sensitivity and
#'       logit specificity Fisher information as \code{size = "fisher"}.
#'       These percentages are then converted to graphical dimensions
#'       using the same RevMan-style normalization applied under
#'       \code{"sampsize_revman"} and \code{"se_revman"}.}
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
#'   Default is \code{TRUE}.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#'  
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
#' Study-specific estimates are shown using the selected symbol (rectangles by default).
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates
#'   \item HSROC curve
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
#' @return
#' No return value. Called for its side effect of producing a plot.
#' @seealso \code{\link{fitRutterGatsonis}}
#' @method plot RutterGatsonis
#' @export
plot.RutterGatsonis <- function(x,
                                symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                                scale=0.02,
                                size=c("fisher","fisher_revman","equal","sampsize","se","sampsize_revman","se_revman"), 
                                specrange=c(0.7,0.995),
                                main="Diagnostic Test Accuracy Meta-Analysis", 
                                HSROC=TRUE,
                                conflevel=0.95,
                                studyCI=FALSE,
                                studylabels=FALSE,
                                ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  size    <- match.arg(size)
  symbol  <- match.arg(symbol)
  warn_unestimable_sroc_points(x$data)
  Lambda  <- x$sdreport2["Lambda","Estimate"]
  beta    <- x$sdreport2["beta","Estimate"]
  roc_points2 <- getROCpoints(Lambda,beta,specrange=specrange)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty = "s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  tryCatch(
  if(size %in% c("fisher","fisher_revman")){
    Y_pw  <- reshapeX_REIT(X=x$data)
    X_pw  <- cbind(Y_pw$sens,Y_pw$spec)
    XT_pw <- t(X_pw)
    nstudy <- nrow(x$data)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    invn[!is.finite(invn)] <- 100000000000 
    A <- diag(invn)
    ####
    Theta  <- x$sdreport2["Theta","Estimate"]
    random <- x$sdreport$par.random
    alpha  <- random[names(random)=="alpha"]
    theta  <- random[names(random)=="theta"]
    lsens  <- (Theta+theta+0.5*(Lambda+alpha))*exp(-0.5*beta)
    lspec  <- -(Theta+theta-0.5*(Lambda+alpha))*exp(0.5*beta)
    p_pw   <- as.vector(rbind(stats::plogis(lsens),
                              stats::plogis(lspec)))
    ####
    var_pw <- p_pw*(1-p_pw)
    B <- diag(var_pw)
    ####
    varA  <- as.numeric(x$Reitsma_recovered["sigma2_A.sens"])
    varB  <- as.numeric(x$Reitsma_recovered["sigma2_B.spec"])
    sAB   <- as.numeric(x$Reitsma_recovered["sigma_AB"])
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
    if(size=="fisher_revman"){
      c_fisher <- 8/mean(pmax(pct$se,pct$sp))
      pct$se   <- (4+c_fisher*pct$se)/10
      pct$sp   <- (4+c_fisher*pct$sp)/10
    }
  }, error = function(e) {
    warning(
      "Fisher-information study sizing could not be calculated. ",
      "Equal study-symbol sizing was used instead. Original error: ",
      conditionMessage(e)) 
    pct <- getWEIGHTS(xdata=x$data,size="equal")})
  if(!size %in% c("fisher","fisher_revman")){
    pct <- getWEIGHTS(xdata=x$data,size=size)
  }
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
  # Add the ROC curve
  if(isTRUE(HSROC)){
    points(roc_points2, type="l", lwd=2,ann=FALSE)###
  # Add the legend 
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "HSROC curve",
                      "Data"), 
           pch = c(NA,NA,symb), 
           lty = c(NA,1,NA), 
           lwd = c(NA,2,NA), 
           col = c(NA,"black","darkgray"))
  }
  invisible(NULL)
}
