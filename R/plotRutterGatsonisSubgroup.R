#' Plot Results from a Rutter and Gatsonis Subgroup Model
#'
#' Produces summary ROC plots for objects of class \code{"RutterGatsonisSubgroup"}
#' obtained from \code{\link{fitRutterGatsonisSubgroup}}. The plot displays
#' study-level estimates of sensitivity and specificity, stratified by subgroup,
#' together with subgroup-specific HSROC (hierarchical summary ROC) curves.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroup"}, as returned by
#'   \code{\link{fitRutterGatsonisSubgroup}}.
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the subgroups. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally. Relative adjustments can be 
#' specified directly, for example, \code{scale = 0.02 * 0.3} produces marker 
#' dimensions equal to 30 percent of the default.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size. Default.}
#'    \item{"fisher"}{Displays parameter-specific percentage study 
#'    contributions to the recovered subgroup-specific pooled logit sensitivity 
#'    and logit specificity. Fisher sizing is available when HSROC accuracy and 
#'    threshold parameters are subgroup-specific.}
#'    \item{"fisher_revman"}{Uses the same parameter-specific percentage Fisher-information
#'      contributions as \code{size = "fisher"}, but applies a fixed-baseline RevMan-style
#'      graphical normalization. Percentage contributions retain their original
#'      interpretation.}
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
#' @param col Vector of colours used for subgroup-specific HSROC curves
#'   and study-level point estimates. If \code{NULL}, colours are generated automatically.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param nudge_legend Numeric horizontal offset for the subgroup legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param connectstudies Whether the point estimates of two subgroups 
#'   within the same study should be connected. Defaults to \code{FALSE}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
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
#'   \item Study-level sensitivity and specificity estimates by subgroup
#'   \item Subgroup-specific HSROC curves
#'   \item A legend identifying the subgroups
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
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitRutterGatsonisSubgroup}}
#' @method plot RutterGatsonisSubgroup
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.RutterGatsonisSubgroup <- function(x, 
                                        symbol=NULL,
                                        scale=0.02,
                                        size=c("equal","fisher","fisher_revman","sampsize","se","sampsize_revman","se_revman"), 
                                        nudge_legend=-0.4,
                                        specrange=c(0.7,0.995),
                                        col=NULL,
                                        main="Diagnostic Test Accuracy Meta-Analysis",
                                        connectstudies=FALSE,
                                        ...){
   if(connectstudies) {
     if(length(unique(x$data$subgroup)) != 2) {
       warning("'connectstudies=TRUE' is only recommended for two-subgroup comparisons." )
     }
   }
   size <- match.arg(size)
   if (size %in% c("fisher","fisher_revman")) {
     if (any(c("accuracy", "threshold") %in% x$constrain)) {
       stop(
         "Fisher study contributions are not currently available when ",
         "HSROC accuracy or threshold effects are constrained to be common ",
         "across subgroups."
       )
     }
   }
   sub  <- x$subgroups
   nsub <- length(sub)
   nstudy <- nrow(x$data)
   warn_unestimable_sroc_points(x$data)
   aes  <- validateSubgroupAesthetics(col=col,symbol=symbol,nsub=nsub)
   col  <- aes$col
   col2 <- aes$col2
   symbols2 <- aes$symbols2
   symb <- vector(mode="integer",length=nsub)
   ##
   oldpar <- par(no.readonly = TRUE)
   on.exit(par(oldpar))
   par(mar = c(5, 4, 4, 10),
       pty="s")   # enlarge right margin
   plot_SESPGRID(main=main)
   # Data points
   if(size %in% c("fisher","fisher_revman")){
     Y_pw  <- reshapeX_REIT(X=x$data)
     X_pw <- matrix(0,nrow=2*nstudy,ncol=2*nsub)
     col_names <- unlist(lapply(sub,function(sg) {c(paste0("mu_A.", sg),paste0("mu_B.", sg))}))
     colnames(X_pw) <- col_names
     se_index <- integer(nstudy)
     sp_index <- integer(nstudy)
     for (i in seq_len(nstudy)) {
       sg <- as.character(x$data$subgroup[i])
       g <- match(sg, sub)
       row_se <- 2 * i - 1
       row_sp <- 2 * i
       se_index[i] <- 2 * g - 1
       sp_index[i] <- 2 * g
       X_pw[row_se, se_index[i]] <- 1
       X_pw[row_sp, sp_index[i]] <- 1
     }
     XT_pw <- t(X_pw)
     nstudy <- nrow(x$data)
     Z <- diag(2*nstudy)
     invn <- 1/Y_pw$n
     invn[!is.finite(invn)] <- 100000000000 
     A <- diag(invn)
     ####
     rn <- rownames(x$sdreport2)
     Lambda <- x$sdreport2[paste0("Lambda_",sub),"Estimate"]
     Theta  <- x$sdreport2[paste0("Theta_",sub),"Estimate"]
     beta   <- x$sdreport2[paste0("beta_",sub),"Estimate"]
     names(Lambda) <- sub("^Lambda\\_","",names(Lambda))
     names(Theta)  <- sub("^Theta\\_","",names(Theta))
     names(beta)   <- sub("^beta\\_","",names(beta))
     Lambda2 <- Lambda[as.character(x$data$subgroup)]
     Theta2  <- Theta[as.character(x$data$subgroup)]
     beta2   <- beta[as.character(x$data$subgroup)]
     random <- x$sdreport$par.random
     alpha  <- random[names(random)=="alpha"]
     theta  <- random[names(random)=="theta"]
     lsens  <- (Theta2+theta+0.5*(Lambda2+alpha))*exp(-0.5*beta2)
     lspec  <- -(Theta2+theta-0.5*(Lambda2+alpha))*exp(0.5*beta2)
     p_pw   <- as.vector(rbind(stats::plogis(lsens),
                               stats::plogis(lspec)))
     ####
     var_pw <- p_pw*(1-p_pw)
     B <- diag(var_pw)
     ####
     G <- matrix(0,ncol=nstudy*2,nrow=nstudy*2)
       for(i in 1:nstudy){
         j       <- 2*i
         sg      <- x$data$subgroup[i]
         varA    <- as.numeric(x$Reitsma_recovered[sg,"sigma2_A.sens"])
         varB    <- as.numeric(x$Reitsma_recovered[sg,"sigma2_B.spec"])
         sAB     <- as.numeric(x$Reitsma_recovered[sg,"sigma_AB"])
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
       se_idx  <- se_index[i]
       sp_idx  <- sp_index[i]
       pct$se[i] <- 100*(weight[se_idx,se_idx]/varb[se_idx,se_idx])
       pct$sp[i] <- 100*(weight[sp_idx,sp_idx]/varb[sp_idx,sp_idx])
     }
     pct$se <- pmax(0,pct$se)
     pct$sp <- pmax(0,pct$sp)
     if(size=="fisher_revman"){
       c_fisher <- 8/mean(pmax(pct$se,pct$sp))
       pct$se   <- (4+c_fisher*pct$se)/10
       pct$sp   <- (4+c_fisher*pct$sp)/10
     }
    } else {
    pct <- getWEIGHTS(xdata=x$data,size=size)
   }
   ###
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
   ###
   lamb <- paste0("Lambda_",sub)
   bet  <- paste0("beta_",sub)
  
   Lambda <- x$sdreport2[lamb,]
   beta   <- x$sdreport2[bet,]
   for(i in seq_len(nsub)){
      roc_points2 <- getROCpoints(Lambda[i,"Estimate"],
                                  beta[i,"Estimate"],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=FALSE,col=col[i])
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
   ###
   legend("right",
          inset = c(nudge_legend, 0),
          legend = sub,
          col    = col,
          pch = symb,
          xpd = TRUE,
          cex = 1.2,
          bty = "n")

   legend("bottomright", 
          bty ="n",
          legend = c(NA,
                     "HSROC curve",
                     "Data"), 
          pch = c(NA,NA,symb[1]), 
          lty = c(NA,1,NA), 
          lwd = c(NA,2,NA), 
          col = c(NA,"black","darkgray"))
   invisible(NULL)
}