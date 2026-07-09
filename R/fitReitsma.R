#' Fit Reitsma Model
#'
#' Fits the Reitsma bivariate random-effects model for diagnostic test accuracy (DTA)
#' meta-analysis using a binomial-normal likelihood via \code{glmmTMB}.
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
#'
#' @return A list of class \code{"Reitsma"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities
#'   \item \code{glmmTMB}: fitted model object.
#'   \item \code{estimates}: parameter estimates with SE.
#'   \item \code{vcov}: variance-covariance matrix.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{LRDOR}: Diagnostic odds ratio and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsma(
#'   data = anticcp,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study
#' )
#' fit$estimates
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
#' @importFrom glmmTMB glmmTMB getME fixef
#' @importFrom stats binomial complete.cases qnorm plogis qlogis vcov cor sd reshape
#' @export
fitReitsma <- function(data,
                       TP, FP, FN, TN,
                       study,
                       constrain=NULL,
                       conflevel=0.95) {
  
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
  XP <- getXP(X=XP)
  
  ### Get initial values
  logit_sens   <- stats::qlogis(pmin(pmax(XP$sens,0.005),0.995))
  logit_spec   <- stats::qlogis(pmin(pmax(XP$spec,0.005),0.995))
  muA_init     <- mean(logit_sens,na.rm=TRUE)
  muB_init     <- mean(logit_spec,na.rm=TRUE)
  sA_init      <- stats::sd(logit_sens,na.rm=TRUE)
  sA_init      <- max(sA_init,1e-05)
  sB_init      <- stats::sd(logit_spec,na.rm=TRUE)
  sB_init      <- max(sB_init,1e-05)
  rAB_init     <- max(min(stats::cor(logit_sens,
                                     logit_spec,
                                     use="pairwise.complete.obs"),0.99),-0.99)
  if(is.na(rAB_init)) rAB_init <- 0
  theta3_init  <- rAB_init/sqrt(1-rAB_init**2)
  
  ###
  start_list   <- list(beta=c(muA_init,muB_init),
                       theta=c(log(sA_init),log(sB_init),theta3_init))
  
  map <- list()
  
  if(!is.null(constrain)){
    if (constrain == "sigma_AB") {
      start_list$theta[3] <- 0
      map$theta <- factor(c(1,2,NA))
    }
    if (constrain == "sigma2_A") {
      start_list$theta[1] <- log(.Machine$double.eps)
      start_list$theta[3] <- 0
      map$theta <- factor(c(NA,1,NA))
    }
    if (constrain == "sigma2_B") {
      start_list$theta[2] <- log(.Machine$double.eps)
      start_list$theta[3] <- 0
      map$theta <- factor(c(1,NA,NA))
    }
    if (constrain == "all") {
      start_list$theta[1] <- log(.Machine$double.eps)
      start_list$theta[2] <- log(.Machine$double.eps)
      start_list$theta[3] <- 0
      map$theta <- factor(c(NA,NA,NA))
    }
  }
  
  ### resphaping the data
  Y    <- reshapeX_REIT(X)
  ### Fitting the Reitsma model
  MA_Y <- glmmTMB::glmmTMB(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec | recordid), 
                           data=Y, family=stats::binomial(link="logit"),
                           start=start_list,
                           map = if(length(map) == 0) NULL else map)
  if (MA_Y$fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  ma_Y      <- summary(MA_Y)
  qq        <- stats::qnorm(1-(1-conflevel)/2)
  ### Sensitivity and Specificity
  sesp           <- as.data.frame(ma_Y$coefficients$cond)
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,(5:8)]
  colnames(sesp) <- c("Estimate","conflevel","CI_Lower","CI_Upper")
  ### SAS variance covariance matrix
  theta      <- glmmTMB::getME(MA_Y,"theta")
  beta_fix   <- glmmTMB::fixef(MA_Y)$cond
  V_full     <- stats::vcov(MA_Y, full = TRUE)
  V_full[is.na(V_full)] <- 0
  esti_V_g   <- getesti_V_g(beta_fix=beta_fix, 
                            theta=theta, 
                            V_full=V_full) 
  esti       <- esti_V_g$esti
  # diagnostic odds ratio, the positive and negative
  # likelihood ratios
  lsens  <- esti[1,"Estimate"]
  lspec  <- esti[2,"Estimate"]
  S      <- ma_Y$vcov$cond
  lrdor  <- getLRDOR(lsens=lsens, lspec=lspec, S=S, conflevel=conflevel)
  # Recover Rutter and Gatsonis estimates
  sigma2_a <- esti[3,"Estimate"]
  sigma2_b <- esti[4,"Estimate"]
  sigma_ab <- esti[5,"Estimate"]
  sigma_a  <- sqrt(sigma2_a)
  sigma_b  <- sqrt(sigma2_b)
  ruga     <- getRUGA(lspec    = lspec,
                      lsens    = lsens,
                      sigma_a  = sigma_a,
                      sigma_b  = sigma_b,
                      sigma_ab = sigma_ab)
  ##
  ret <- list(data      = XP,
              glmmTMB   = MA_Y,
              estimates = esti,
              vcov      = esti_V_g$V_g,
              sensspec  = sesp,
              LRDOR     = lrdor,
              RutterGatsonis_recovered = ruga,
              constrain = constrain)
  class(ret) <- c("Reitsma","Cochrane")
  return(ret)
}