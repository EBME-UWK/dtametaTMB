#' Fit the Reitsma LCA Model
#'
#' Fits the Reitsma latent class model for diagnostic test accuracy (DTA),
#' allowing for an imperfect reference standard under conditional independence.
#'
#' @param data A data.frame containing study-level data.
#' @param y11 Number testing positive with both index and reference tests (column name).
#' @param y10 Number testing positive with index test but negative with reference test (column name).
#' @param y01 Number testing negative with index test but positive with reference test (column name).
#' @param y00 Number testing negative with both index and reference tests (column name).
#' @param study Study identifier (column name). 
#' @param constrain Optional character string specifying a simplified
#'   covariance structure for the index test of the Reitsma LCA model.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       The standard unconstrained Reitsma model is fitted.
#'     }
#'     \item{\code{"sigma_AB.index"}}{
#'       The covariance between logit-sensitivity and
#'       logit-specificity random effects is fixed at zero.
#'       Random effects remain independent.
#'     }
#'     \item{\code{"sigma2_A.index"}}{
#'       The between-study variance of logit-sensitivity is fixed at zero.
#'       This also implies a zero covariance.
#'     }
#'     \item{\code{"sigma2_B.index"}}{
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
#' @return A list of class \code{"ReitsmaLCA"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities.
#'   \item \code{fit}: Optimization result from \code{nlminb}.
#'   \item \code{sdreport}: TMB standard report.
#'   \item \code{sdreport2}: Summary of reported parameters.
#'   \item \code{vcov}: variance-covariance matrix of reported parameters.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{prevref}: Estimated (average) prevalence and reference standard sensitivity/specificity with confidence intervals.
#'   \item \code{LRDOR}: Diagnostic odds ratio and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsmaLCA(
#'   data = anticcp,
#'   y11 = TP,
#'   y10 = FP,
#'   y01 = FN,
#'   y00 = TN,
#'   study = study
#' )
#' fit
#' summary(fit)
#' 
#'
#' @references
#' Liu, Y., Chen, Y., & Chu, H. (2015). 
#' A unification of models for meta-analysis of diagnostic accuracy studies without a gold standard. 
#' \emph{Biometrics}, 71(2), 538-547.
#' \doi{10.1111/biom.12264}
#' 
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
#' @importFrom stats complete.cases qnorm plogis qlogis nlminb var
#' @export

fitReitsmaLCA <- function(data,
                          y11,
                          y10,
                          y01,
                          y00,
                          study,
                          conflevel=0.95,
                          constrain=NULL,
                          verbose=FALSE){
  
  allowed_constraints <- c(
    "sigma_AB.index",
    "sigma2_A.index",
    "sigma2_B.index",
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
  
  y11_col <- deparse(substitute(y11))
  y10_col <- deparse(substitute(y10))
  y01_col <- deparse(substitute(y01))
  y00_col <- deparse(substitute(y00))
  study_col <- deparse(substitute(study))
  
  dat <- data.frame(
    study = data[[study_col]],
    y11 = data[[y11_col]],
    y10 = data[[y10_col]],
    y01 = data[[y01_col]],
    y00 = data[[y00_col]]
  )
  
  X <- check2_data(dat=dat,conflevel=conflevel)
  X$n     <- with(X, y11+y10+y01+y00)
  n_study <- nrow(X)
  
  ### Get initial values
  init <- tryCatch(fitReitsmaTMB(data=X,
                                 TP=y11,FP=y10,FN=y01,TN=y00,study=study,
                                 constrain=NULL,
                                 conflevel=conflevel)$sdreport,
                   error = function(e) NULL)
  muA_init <- NA_real_
  muB_init <- NA_real_
  lsA_init <- NA_real_
  lsB_init <- NA_real_
  thetaAB_init <- NA_real_
  mup_init <- NA_real_
  s2p_init <- NA_real_
  muA_ref_init <- NA_real_
  muB_ref_init <- NA_real_
  if(!is.null(init)){
    muA_init     <- init$par.fixed["mu_A"]
    muB_init     <- init$par.fixed["mu_B"]
    lsA_init     <- init$par.fixed["log_sigma_A"]
    lsB_init     <- init$par.fixed["log_sigma_B"]
    thetaAB_init <- init$par.fixed["theta_AB"]
  }
  if (!is.finite(muA_init)) {muA_init <- 0}
  if (!is.finite(muB_init)) {muB_init <- 0}
  if (!is.finite(lsA_init)) {lsA_init <- log(0.5)}
  if (!is.finite(lsB_init)) {lsB_init <- log(0.5)}
  if (!is.finite(thetaAB_init)) {thetaAB_init <- 0}
  
  muA_ref_init <- stats::qlogis(mean(c(stats::plogis(muA_init),0.99)))
  muB_ref_init <- stats::qlogis(mean(c(stats::plogis(muB_init),0.99)))
  
  lprev    <- stats::qlogis(with(X,(y11+y01+0.5)/(y11+y10+y01+y00+1)))
  mup_init <- mean(lprev)
  if (!is.finite(mup_init)) {mup_init <- 0}
  s2p_init <- stats::var(lprev)
  if (!is.finite(s2p_init) || s2p_init <= 0) {s2p_init <- 0.25}
  s2p_init <- max(s2p_init,1e-10)
  
  parameters <- list(
    mu_A_index = muA_init,
    mu_B_index = muB_init,
    mu_A_ref   = muA_ref_init,
    mu_B_ref   = muB_ref_init,
    mu_prev    = mup_init,
    
    log_sigma_prev     = 0.5*log(s2p_init),
    log_sigma_A_index  = lsA_init,
    log_sigma_B_index  = lsB_init,
    theta_AB_index     = thetaAB_init,
    
    prevu = rep(0,n_study),
    sensu = rep(0,n_study),
    specu = rep(0,n_study)
  )
  
  dat2 <- list(
    y11 = X$y11,
    y10 = X$y10,
    y01 = X$y01,
    y00 = X$y00
  )
  
  #################
  ### Constrain ###
  #################
  
  map <- list()
  
  if(!is.null(constrain)){
    if (constrain == "sigma_AB.index") {
      parameters$theta_AB_index <- 0
      map$theta_AB_index        <- factor(NA)
    }
    if (constrain == "sigma2_A.index") {
      parameters$log_sigma_A_index <- log(.Machine$double.eps)
      parameters$theta_AB_index    <- 0
      map$log_sigma_A_index        <- factor(NA)
      map$theta_AB_index           <- factor(NA)
    }
    if (constrain == "sigma2_B.index") {
      parameters$log_sigma_B_index <- log(.Machine$double.eps)
      parameters$theta_AB_index    <- 0
      map$log_sigma_B_index        <- factor(NA)
      map$theta_AB_index           <- factor(NA)
    }
    if (constrain == "all") {
      parameters$log_sigma_A_index <- log(.Machine$double.eps)
      parameters$log_sigma_B_index <- log(.Machine$double.eps)
      parameters$theta_AB_index    <- 0
      map$log_sigma_A_index        <- factor(NA)
      map$log_sigma_B_index        <- factor(NA)
      map$theta_AB_index           <- factor(NA)
    }
  }

  dat2$model <- "ReitsmaLCA"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("prevu", "sensu", "specu"),
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
  
  # Reports
  rep  <- TMB::sdreport(obj)
  rep2 <- summary(rep, select = "report")
  
  rownames(rep2) <- gsub("_ref", ".ref",
                         gsub("_index", ".index",
                              rownames(rep2)))
  
  # Recover Rutter and Gatsonis Parameters
  indexRUGA <- getRUGA(lsens=rep2["mu_A.index","Estimate"],
                       lspec=rep2["mu_B.index","Estimate"],
                       sigma_a=sqrt(rep2["sigma2_A.index","Estimate"]),
                       sigma_b=sqrt(rep2["sigma2_B.index","Estimate"]),
                       sigma_ab=rep2["sigma_AB.index","Estimate"])

  # Get variance covariance matrix
  vcov <- rep$cov
  colnames(vcov) <- rownames(vcov) <- rownames(rep2)
  
  # Get empirical Bayes estimates
  u_est <- rep$par.random
  u_var <- rep$diag.cov.random
  u     <- data.frame(raneff=names(u_est),
                      u=u_est,
                      u_var=u_var)
  
  lsens_study <- rep2["mu_A.index","Estimate"] + u[u$raneff=="sensu","u"]
  lspec_study <- rep2["mu_B.index","Estimate"] + u[u$raneff=="specu","u"]
  
  X$sens_eb  <- stats::plogis(lsens_study)
  X$spec_eb  <- stats::plogis(lspec_study)
  
  X$lsens_eb_var <- u[u$raneff=="sensu","u_var"]
  X$lspec_eb_var <- u[u$raneff=="specu","u_var"]
  
  # How to get sensitivities
  qq   <- stats::qnorm(1-(1-conflevel)/2)
  ### Sensitivity and Specificity
  sesp           <- as.data.frame(rep2[c("mu_prev",
                                         "mu_A.index",
                                         "mu_B.index",
                                         "mu_A.ref",
                                         "mu_B.ref"),])
  sesp$type      <- c("Prev","Sens","Spec","Sens","Spec")
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,c("type","Orig","conflevel","CI_Lower","CI_Upper")]
  colnames(sesp) <- c("type","Estimate","conflevel","CI_Lower","CI_Upper")
  sesp1          <- sesp[c("mu_A.index","mu_B.index"),]
  prevref        <- sesp[c("mu_prev","mu_A.ref","mu_B.ref"),]
  ####
  # likelihood ratios
  lsens  <- rep2["mu_A.index","Estimate"]
  lspec  <- rep2["mu_B.index","Estimate"]
  S      <- vcov[c("mu_A.index","mu_B.index"),
                 c("mu_A.index","mu_B.index")]
  lrdor  <- getLRDOR(lsens=lsens, lspec=lspec, S=S, conflevel=conflevel)
  lrdor$test <- "index"
  lrdor  <- lrdor[,c("test","Estimate","conflevel","CI_Lower","CI_Upper")]
  # Result object
  res <- list(
    data         = X,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    vcov         = vcov,
    sensspec     = sesp1,
    prevref      = prevref,
    LRDOR        = lrdor,
    RutterGatsonis_recovered = indexRUGA,
    constrain    = constrain
  )
  
  # Assign class
  class(res) <- c("ReitsmaLCA","Reitsma","DTAmodel")
  
  return(res)
}