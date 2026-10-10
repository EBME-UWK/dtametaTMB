#' Fit the Reitsma TMB Model
#'
#' Fits the Reitsma bivariate random-effects model for diagnostic test accuracy (DTA)
#' meta-analysis using a binomial-normal likelihood via a custom \code{TMB} template,
#' which may offer improved numerical robustness in difficult or near-boundary cases.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name). 
#' @param constrain Optional character string specifying a simplified
#'   covariance structure for the Reitsma model.
#'
#'   This can be useful for sparse data, small meta-analyses, or for
#'   reproducing simplified bivariate models described in the Cochrane
#'   Handbook for Diagnostic Test Accuracy Reviews.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       The standard unconstrained Reitsma model is fitted.
#'     }
#'     \item{\code{"sigma_AB"}}{
#'       The covariance between logit-sensitivity and
#'       logit-specificity random effects is fixed at zero.
#'       Random effects remain independent.
#'     }
#'     \item{\code{"sigma2_A"}}{
#'       The between-study variance of logit-sensitivity is fixed at zero.
#'       This also implies a zero covariance.
#'     }
#'     \item{\code{"sigma2_B"}}{
#'       The between-study variance of logit-specificity is fixed at zero.
#'       This also implies a zero covariance.
#'     }
#'     \item{\code{"all"}}{
#'       All random-effects variance and covariance parameters are fixed
#'       at zero, resulting in a fixed-effects model.
#'     }
#'   }
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return A list of class \code{"ReitsmaTMB"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities
#'   \item \code{fit}: Optimization result from \code{nlminb}.
#'   \item \code{sdreport}: TMB standard report.
#'   \item \code{sdreport2}: Summary of reported parameters.
#'   \item \code{vcov}: variance-covariance matrix of reported parameters.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{LRDOR}: Diagnostic odds ratio and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsmaTMB(
#'   data = anticcp,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study
#' )
#' summary(fit)
#' 
#'
#' @references
#' Reitsma, J. B., et al. (2005). 
#' Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews.
#' \emph{Journal of Clinical Epidemiology}, 58(10), 982–990.
#' \doi{10.1016/j.jclinepi.2005.02.022}
#'
#' Rutter, C. M., & Gatsonis, C. A. (2001). 
#' A hierarchical regression approach to meta-analysis of diagnostic test accuracy evaluations.
#' \emph{Statistics in Medicine}, 20(19), 2865–2884.
#' \doi{10.1002/sim.942}
#' 
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats binomial complete.cases qnorm plogis qlogis vcov cor sd reshape
#' @export
fitReitsmaTMB <- function(data,
                          TP, FP, FN, TN,
                          study,
                          constrain=NULL,
                          conflevel=0.95,
                          verbose=FALSE) {
  
  allowed_constraints <- c(
    "sigma_AB",
    "sigma2_A",
    "sigma2_B",
    "all"
  )
  
  if (!is.null(constrain)) {
    if (!is.character(constrain) ||
        length(constrain) != 1 ||
        !constrain %in% allowed_constraints) {
      stop(
        "'constrain' must be one of: ",
        paste(shQuote(allowed_constraints), collapse = ", "),
        " or NULL."
      )
    }
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  
  TP_col <- deparse(substitute(TP))
  FP_col <- deparse(substitute(FP))
  FN_col <- deparse(substitute(FN))
  TN_col <- deparse(substitute(TN))
  study_col <- deparse(substitute(study))
  
  dat <- data.frame(
    study = data[[study_col]],
    TP = data[[TP_col]],
    TN = data[[TN_col]],
    FP = data[[FP_col]],
    FN = data[[FN_col]]
  )
  
  X <- XP <- check_data(dat=dat,
                        conflevel=conflevel)
  n_study <- nrow(X)
  XP <- getXP(X=XP)
  
  ### Get initial values
  init <- initReitsma(X=X)
  ###
  parameters <- list(
    mu_A        = init$muA_init,
    mu_B        = init$muB_init,
    log_sigma_A = log(init$sA_init),
    log_sigma_B = log(init$sB_init),
    theta_AB    = atanh(init$rAB_init),
    sensu = rep(0,n_study),
    specu = rep(0,n_study)
  )
  
  map <- list()
  
  
  if(!is.null(constrain)){
    if (constrain == "sigma_AB") {
      parameters$theta_AB <- 0
      map$theta_AB        <- factor(NA)
    }
    if (constrain == "sigma2_A") {
      parameters$log_sigma_A <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_A        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
    if (constrain == "sigma2_B") {
      parameters$log_sigma_B <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_B        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
    if (constrain == "all") {
      parameters$log_sigma_A <- log(.Machine$double.eps)
      parameters$log_sigma_B <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_A        <- factor(NA)
      map$log_sigma_B        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
  }
  
  dat2 <- list(
    TP = X$TP,
    FP = X$FP,
    FN = X$FN,
    TN = X$TN
  )
  
  dat2$model <- "Reitsma"

  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("sensu","specu"),
                        silent = !verbose,
                        DLL = "dtametaTMB_TMBExports")

  # Optimization
  fit <- stats::nlminb(obj$par,
                       obj$fn,
                       obj$gr)
  # Convergence warning.
  if (fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  
  if (!sdreport(obj)$pdHess) {
    warning(
      "The Hessian matrix is not positive definite. ",
      "Parameter estimates and standard errors may be unreliable. ",
      "Consider checking model specification, data quality, or potential identifiability issues."
    )
  }
  
  # Reports
  rep  <- TMB::sdreport(obj)
  rep2 <- summary(rep, select = "report")
  
  rownames(rep2)[1:4] <- c("mu_A.sens",
                           "mu_B.spec",
                           "sigma2_A.sens",
                           "sigma2_B.spec")
  
  # Recover Rutter and Gatsonis Parameters
  RUGA <- getRUGA(lsens=rep2["mu_A.sens","Estimate"],
                  lspec=rep2["mu_B.spec","Estimate"],
                  sigma_a=sqrt(rep2["sigma2_A.sens","Estimate"]),
                  sigma_b=sqrt(rep2["sigma2_B.spec","Estimate"]),
                  sigma_ab=rep2["sigma_AB","Estimate"])
  # Get variance covariance matrix
  vcov <- rep$cov
  colnames(vcov) <- rownames(vcov) <- rownames(rep2)
  #####
  qq        <- stats::qnorm(1-(1-conflevel)/2)
  ### Sensitivity and Specificity
  sesp           <- as.data.frame(rep2[c("mu_A.sens","mu_B.spec"),],
                                  row.names=c("sens","spec"))
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,c("Orig","conflevel","CI_Lower","CI_Upper")]
  colnames(sesp) <- c("Estimate","conflevel","CI_Lower","CI_Upper")
  # diagnostic odds ratio, the positive and negative
  # likelihood ratios
  lsens  <- rep2["mu_A.sens","Estimate"]
  lspec  <- rep2["mu_B.spec","Estimate"]
  S      <- vcov[c("mu_A.sens","mu_B.spec"),c("mu_A.sens","mu_B.spec")]
  lrdor  <- getLRDOR(lsens=lsens, lspec=lspec, S=S, conflevel=conflevel)
  ##
  res <- list(
    data         = XP,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    vcov         = vcov,
    sensspec     = sesp,
    LRDOR        = lrdor,
    RutterGatsonis_recovered = RUGA,
    constrain    = constrain
  )

  # Assign class
  class(res) <- c("ReitsmaTMB","Reitsma","DTAmodel")

  return(res)
}


#' Print Method for ReitsmaTMB Objects
#'
#' Displays a concise summary of a fitted ReitsmaTMB diagnostic test
#' accuracy model, including number of studies, convergence status
#' and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"ReitsmaTMB"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.ReitsmaTMB}}
#' @return
#' Invisibly returns the input object.
#' @method print ReitsmaTMB
#' @export
print.ReitsmaTMB <- function(x, ...) {
  
  cat("\n", "Reitsma TMB Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma TMB Model")), "\n\n", sep = "")
  
  
  n_study   <- nrow(x$data)
  converged <- x$fit$convergence == 0
  pdHess    <- x$sdreport$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies :", n_study, "\n")
  cat("Optimizer         :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian           :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|        :", max(abs(x$sdreport$gradient.fixed)), "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n")
  
  cat("\n")
  
  est <- x$sdreport2[,"Estimate"]
  
  cat("mu_A (logitsens)  :", round(est["mu_A.sens"], 3), "\n")
  cat("mu_B (logitspec)  :", round(est["mu_B.spec"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.ReitsmaTMB <- logLik.DTAmodel

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaTMB <- coef.DTAmodel

#' Summary Method for ReitsmaTMB Objects
#'
#' Extracts key results from an object of class \code{"ReitsmaTMB"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"ReitsmaTMB"} as returned by
#'   \code{\link{fitReitsmaTMB}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class \code{"summary.ReitsmaTMB"} with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#' }
#'
#' @seealso \code{\link{fitReitsmaTMB}}
#'
#'
#' @export
summary.ReitsmaTMB <- function(object, ...) {
  ret <- list(estimates = object$sdreport2,
              sensspec = object$sensspec,
              RutterGatsonis_recovered = object$RutterGatsonis_recovered)
  class(ret) <- c("summary.ReitsmaTMB","summary.Reitsma","summary.DTAmodel")
  return(ret)
}



#' Plot Results from a ReitsmaTMB Model
#'
#' Produces a summary ROC plot for objects of class \code{"ReitsmaTMB"}
#' obtained from \code{\link{fitReitsmaTMB}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. Optionally, the
#' HSROC (hierarchical summary ROC) curve can be overlaid.
#'
#' @param x An object of class \code{"ReitsmaTMB"}, as returned by
#'   \code{\link{fitReitsmaTMB}}.
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
#'    \item{"fisher"}{Size proportional to a decomposition of Fisher's Information matrix. Default.}
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
#' @param HSROC If \code{TRUE}, the HSROC curve is added to the plot.
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
#' @param studyCI Whether study-level confidence
#'   intervals should be displayed for the observed sensitivity and
#'   specificity estimates. Defaults to \code{FALSE}.
#'
#' @param studylabels Whether study identifiers should
#'   be displayed next to the study-level estimates. Defaults to
#'   \code{FALSE}.
#'
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
#' @return
#' No return value. Called for its side effect of producing a plot.
#' @seealso \code{\link{fitReitsmaTMB}}
#' @importFrom stats qlogis plogis predict qf
#' @method plot ReitsmaTMB
#' @export
plot.ReitsmaTMB <- function(x,
                            symbol=c("rectangle","ellipse","diamond","triangle","cross","plus","star"),
                            scale=0.02,
                            size=c("fisher","fisher_revman","equal","sampsize","se","sampsize_revman","se_revman"), 
                            main="Diagnostic Test Accuracy Meta-Analysis",
                            HSROC=FALSE, 
                            specrange=c(0.7,0.995),
                            conflevel=0.95,
                            predlevel=0.95, 
                            studyCI=FALSE,
                            studylabels=FALSE,
                            ...) {
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
  muA     <- x$sdreport2["mu_A.sens","Estimate"]
  muB     <- x$sdreport2["mu_B.spec","Estimate"]
  seB     <- x$sdreport2["mu_B.spec","Std. Error"]
  seA     <- x$sdreport2["mu_A.sens","Std. Error"]
  covAB   <- x$vcov["mu_A.sens","mu_B.spec"]
  varA    <- x$sdreport2["sigma2_A.sens","Estimate"]
  varB    <- x$sdreport2["sigma2_B.spec","Estimate"]
  sAB     <- x$sdreport2["sigma_AB","Estimate"]
  region <- getConfPredRegion(muA=muA,muB=muB,
                              seA=seA,seB=seB,covAB=covAB, # conf
                              varA=varA,varB=varB,sAB=sAB, # pred
                              nstudy=nstudy,
                              conflevel=conflevel,
                              predlevel=predlevel)
  # Calculations for percentage weights
  if(size %in% c("fisher","fisher_revman")){
  pct <- tryCatch(
  {   
    Y_pw  <- reshapeX_REIT(X=x$data)
    X_pw  <- cbind(Y_pw$sens, Y_pw$spec)
    XT_pw <- t(X_pw)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    invn[!is.finite(invn)] <- 100000000000 
    A <- diag(invn)
    p_pwTMB  <- x$sdreport$par.random
    p_pwTMB2 <- data.frame(sens=stats::plogis(muA+p_pwTMB[names(p_pwTMB)=="sensu"]),
                           spec=stats::plogis(muB+p_pwTMB[names(p_pwTMB)=="specu"]))
    p_pw <- as.vector(rbind(p_pwTMB2$sens,p_pwTMB2$spec))
    ###
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
    if(size=="fisher_revman"){
      c_fisher <- 8/mean(pmax(pct$se,pct$sp))
      pct$se   <- (4+c_fisher*pct$se)/10
      pct$sp   <- (4+c_fisher*pct$sp)/10
    }
    if (!all(is.finite(as.matrix(pct)))) stop("non-finite Fisher weights")
    pct
  }, error = function(e) {
    warning(
      "Fisher-information study sizing could not be calculated. ",
      "Equal study-symbol sizing was used instead. Original error: ",
      conditionMessage(e)) 
    getWEIGHTS(xdata=x$data,size="equal")})
  } else {
    pct <- getWEIGHTS(xdata=x$data,size=size)
  }
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
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
  if(isTRUE(HSROC)){
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




#' @rdname as_revman.dtametaTMB
#' @export
as_revman.ReitsmaTMB <- function(x, ...) {
  
  ## extract these from x
  mu_se   <- x$sdreport2["mu_A.sens","Estimate"]
  mu_sp   <- x$sdreport2["mu_B.spec","Estimate"]
  var_se  <- x$sdreport2["sigma2_A.sens","Estimate"]
  var_sp  <- x$sdreport2["sigma2_B.spec","Estimate"]
  cov_ss  <- x$sdreport2["sigma_AB","Estimate"]
  cor_ss  <- x$sdreport2["rho_AB","Estimate"]
  
  Lambda  <- x$RutterGatsonis_recovered$Lambda
  Theta   <- x$RutterGatsonis_recovered$Theta
  beta    <- x$RutterGatsonis_recovered$beta
  varA    <- x$RutterGatsonis_recovered$sigma2_alpha
  varT    <- x$RutterGatsonis_recovered$sigma2_theta
  
  seelse  <- x$sdreport2["mu_A.sens","Std. Error"]
  seelsp  <- x$sdreport2["mu_B.spec","Std. Error"]
  coves   <- x$vcov["mu_A.sens","mu_B.spec"]
  nstudy  <- nrow(x$data)
  
  ret <- data.frame(
    Externally_Calculated_Parameters=
      c(rep("HSROC model parameters",5),
        rep("Bivariate model parameters",6),
        rep("Confidence and prediction regions",4)),
    Parameter = c(
      "Lambda",
      "Theta",
      "beta",
      "Var(accuracy)",
      "Var(threshold)",
      "E(logitSe)",
      "E(logitSp)",
      "Var(logitSe)",
      "Var(logitSp)",
      "Cov(logits)",
      "Corr(logits)",
      "SE(E(logitSe))",
      "SE(E(logitSp))",
      "Cov(Es)",
      "Studies"
    ),
    Estimate = c(Lambda,Theta,beta,varA,varT,
                 mu_se,mu_sp,var_se,var_sp,cov_ss,cor_ss,
                 seelse,seelsp,coves,nstudy), 
    
    row.names = NULL
  )
  return(ret)
}
