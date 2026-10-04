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
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the subgroups. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"fisher"}{Study-symbol dimensions represent percentage contributions 
#'      to the estimated logit sensitivity and logit specificity of the study’s own subgroup.
#'      Weights sum to 100 separately within each subgroup and outcome dimension.}
#'    \item{"equal"}{All studies shown with equal size.}
#'    \item{"sampsize"}{
#'      Horizontal and vertical symbol dimensions are proportional to the
#'      relative numbers of non-diseased and diseased participants,
#'      respectively.}
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
#'
#' @seealso \code{\link{fitReitsmaSubgroup}}
#' @method plot ReitsmaSubgroup
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.ReitsmaSubgroup <- function(x,
                                 symbol=NULL,
                                 scale=0.02, 
                                 size=c("fisher","equal","sampsize","se","sampsize_revman","se_revman"), 
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
  warn_unestimable_sroc_points(x$data)
  if(!all(subs == make.names(sub))){
    stop("Object is corrupted. Please don't change object after running fitReitsmaSubgroup().")
  }
  if(!all(x$data$subgroup_safe == make.names(x$data$subgroup))){
    stop("Object is corrupted. Please don't change object after running fitReitsmaSubgroup().")
  }
  if(is.null(col)) col <- grDevices::rainbow(n=nsub)
  col2 <- grDevices::adjustcolor(col,alpha.f=0.6)
  # symbols
  if(is.null(symbol)) {
    available_symbols <- c("rectangle","ellipse","diamond","triangle","cross","plus","star")
    symbols2 <- rep(available_symbols,length.out=nsub)
  } else {
    symbols2 <- symbol
  }
  symb <- vector(mode="integer",length=nsub)
  # Calculations for percentage weights
  if(size=="fisher"){
    X <- x$data
    nstudy  <- nrow(X)
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
    if(is.null(x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=2*nsub)
      col_names <- unlist(lapply(subs,function(sg) {c(paste0("mu_A.", sg),paste0("mu_B.", sg))}))
      colnames(X_pw) <- col_names
      se_index <- integer(nstudy)
      sp_index <- integer(nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup_safe[i])
        g <- match(sg, subs)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        se_index[i] <- 2 * g - 1
        sp_index[i] <- 2 * g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("sens" %in% x$sensspec_constrain & !("spec" %in% x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=nsub+1)
      colnames(X_pw) <- c("mu_A.common",paste0("mu_B.", subs))
      se_index <- rep(1,nstudy)
      sp_index <- integer(nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup_safe[i])
        g  <- match(sg, subs)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        sp_index[i] <- 1 + g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("spec" %in% x$sensspec_constrain & !("sens" %in% x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=nsub+1)
      colnames(X_pw) <- c(paste0("mu_A.", subs),"mu_B.common")
      se_index <- integer(nstudy)
      sp_index <- rep(nsub+1,nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup_safe[i])
        g  <- match(sg, subs)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        se_index[i] <- g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("sens" %in% x$sensspec_constrain & "spec" %in% x$sensspec_constrain){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=2)
      colnames(X_pw) <- c("mu_A.common","mu_B.common")
      se_index <- rep(1,nstudy)
      sp_index <- rep(2,nstudy)
      for (i in seq_len(nstudy)) {
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        X_pw[row_se, 1] <- 1
        X_pw[row_sp, 2] <- 1
      } 
    }  
    XT_pw <- t(X_pw)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    invn[!is.finite(invn)] <- 100000000000 
    A <- diag(invn)
    p_pw <- stats::predict(x$glmmTMB_mu, type="response")
    var_pw <- p_pw*(1-p_pw)
    B <- diag(var_pw)
    G <- matrix(0,ncol=nstudy*2,nrow=nstudy*2)
      for(i in 1:nstudy){
        j       <- 2*i
        sg      <- x$data$subgroup_safe[i]
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
        varA    <- x$estimates_mu[s2_A.sg,"Estimate"]
        varB    <- x$estimates_mu[s2_B.sg,"Estimate"]
        sAB     <- x$estimates_mu[s_AB.sg,"Estimate"]
        G[(j-1):j,(j-1):j] <- matrix(c(varA,sAB,sAB,varB),2,2)
      }
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
      se_idx  <- se_index[i]
      sp_idx  <- sp_index[i]
      pct$se[i] <- 100*(weight[se_idx,se_idx]/varb[se_idx,se_idx])
      pct$sp[i] <- 100*(weight[sp_idx,sp_idx]/varb[sp_idx,sp_idx])
    }
    pct$se <- pmax(0,pct$se)
    pct$sp <- pmax(0,pct$sp)
  } else {
  pct <- getWEIGHTS(xdata=x$data,size=size)
  }
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  for (i in seq_along(sub)){
    pointsXY(x=1-x$data$spec[x$data$subgroup==sub[i]], 
             y=x$data$sens[x$data$subgroup==sub[i]], 
             symbol = symbols2[i], 
             scale = scale*0.5,
             cex.x = pct$sp[x$data$subgroup==sub[i]],
             cex.y = pct$se[x$data$subgroup==sub[i]],
             col=col2[i])
    symb[i] <- switch(symbols2[i], rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
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
    if("sens" %in% x$sensspec_constrain){
      mu_A.sg <- paste0("mu_A.",subs[1])
      } else {
      mu_A.sg <- paste0("mu_A.",sg)
    }
    if("spec" %in% x$sensspec_constrain){
      mu_B.sg <- paste0("mu_B.",subs[1])
      } else {
      mu_B.sg <- paste0("mu_B.",sg) 
    }
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
                      "Data"), 
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
