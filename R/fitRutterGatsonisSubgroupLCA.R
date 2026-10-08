#' Fit the Rutter and Gatsonis Subgroup LCA Model
#'
#' The function fits an HSROC model with a single categorical study-level
#' covariate (subgroup), allowing for an imperfect reference standard assuming
#' conditional independence.
#'
#' @param data A data.frame containing study-level data.
#' @param y11 Number testing positive with both index and reference tests (column name).
#' @param y10 Number testing positive with index test but negative with reference test (column name).
#' @param y01 Number testing negative with index test but positive with reference test (column name).
#' @param y00 Number testing negative with both index and reference tests (column name).
#' @param study Study identifier (column name).
#' @param subgroup A single categorical study-level subgroup variable (column name).
#' @param constrain Optional character vector specifying model parameters
#'   that should be constrained during estimation.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{"sigma2_alpha"}{
#'       Fix the between-study variance of the HSROC accuracy random effect
#'       at zero.
#'     }
#'     \item{"sigma2_theta"}{
#'       Fix the between-study variance of the HSROC threshold random effect
#'       at zero.
#'     }
#'     \item{"accuracy"}{
#'       No subgroup effect on accuracy.
#'     }
#'     \item{"threshold"}{
#'       No subgroup effect on threshold.
#'     }
#'     \item{"shape"}{
#'       No subgroup effect on shape.
#'     }
#'     \item{"shape_zero"}{
#'       All shape parameters (beta's) are fixed at zero.
#'     }
#'   }
#'
#'   Multiple constraints may be specified simultaneously, e.g.
#'   \code{constrain = c("sigma2_alpha", "shape")}.
#'   
#'   \code{shape_zero} overrides \code{shape}.
#'
#'   If \code{NULL} (default), the unconstrained HSROC model is fitted.
#' @param spec Optional specificity value at which sensitivity is estimated.
#'    If \code{NULL}, a specificity of 0.8 is used as a proxy.
#' @param prev_variances Whether the between-study random-effects variance for 
#'   logit prevalence should be assumed to be \code{"common"} (default) or \code{"unequal"} 
#'   across subgroups. If \code{"common"}, a single between-study variance is estimated and shared across all subgroups. 
#'   If \code{"unequal"}, subgroup-specific variances are estimated.
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return
#' An object of class \code{"RutterGatsonisSubgroupLCA"} containing:
#' \describe{
#'   \item{data}{Processed input data with derived quantities.}
#'   \item{fit}{Optimization result from \code{nlminb}.}
#'   \item{sdreport}{TMB standard report.}
#'   \item{sdreport2}{Summary of reported subgroup-specific parameters.}
#'   \item{sensspec}{Estimated subgroup-specific index test sensitivities at the given
#'   specificity value(s), with confidence intervals.}
#'   \item{prevref}{Estimated (average) prevalence and reference standard sensitivity/specificity with confidence intervals.}
#'   \item{constrain}{Constraints on parameters applied during model fitting.}
#'   \item{subgroups}{The subgroup levels used in the model fit.}
#'   \item{prev_variances}{Prevalence variance structure used in the fitted model.}
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats complete.cases median nlminb qnorm plogis
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
#' @examples
#' data("schuetz")
#' fit <- fitRutterGatsonisSubgroupLCA(
#'   data = schuetz,
#'   y11 = TP,
#'   y10 = FP,
#'   y01 = FN,
#'   y00 = TN,
#'   study = study,
#'   subgroup = test
#' )
#' fit
#' summary(fit)
#'
#' @note Requires a compiled TMB model named \code{"RutterGatsonisSubgroupLCA"}.
#' @export
fitRutterGatsonisSubgroupLCA <- function(data,
                                         y11, 
                                         y10, 
                                         y01, 
                                         y00,
                                         study,
                                         subgroup,
                                         constrain=NULL,
                                         prev_variances=c("common","unequal"),
                                         spec=NULL,
                                         conflevel=0.95,
                                         verbose=FALSE){
  
  prev_variances  <- match.arg(prev_variances)
  
  allowed_constraints <- c(
    "sigma2_alpha",
    "sigma2_theta",
    "accuracy",
    "threshold",
    "shape",
    "shape_zero"
  )
  
  if (!is.null(constrain)) {
    
    if (!is.character(constrain)) {
      stop("'constrain' must be a character vector or NULL.")
    }
    
    invalid_constraints <- setdiff(
      constrain,
      allowed_constraints
    )
    
    if (length(invalid_constraints) > 0) {
      stop(
        "Unknown constraint(s): ",
        paste(invalid_constraints,
              collapse = ", ")
      )
    }
    
  }
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  
  y11_col       <- deparse(substitute(y11))
  y10_col       <- deparse(substitute(y10))
  y01_col       <- deparse(substitute(y01))
  y00_col       <- deparse(substitute(y00))
  study_col     <- deparse(substitute(study))
  subgroup_col  <- deparse(substitute(subgroup))
  
  X <- data.frame(study = data[[study_col]],
                  y11 = data[[y11_col]],
                  y10 = data[[y10_col]],
                  y01 = data[[y01_col]],
                  y00 = data[[y00_col]],
                  subgroup = data[[subgroup_col]])
  
  X <- check2_data(dat=X,conflevel=conflevel)
  
  X$subgroup <- droplevels(factor(X$subgroup))
  lsub       <- levels(X$subgroup)
  n_study    <- nrow(X)
  llsub      <- length(lsub)
  X$n        <- with(X, y11+y10+y01+y00)
  
  # Get starting values
  init <- tryCatch(fitRutterGatsonisLCA(data=X,
                                        y11=y11,
                                        y10=y10,
                                        y01=y01,
                                        y00=y00,
                                        study=study,
                                        conflevel=conflevel,
                                        constrain=NULL)$sdreport,
                   error = function(e) NULL)
  
  Lambda_init <- NA_real_
  Theta_init <- NA_real_
  beta_init <- NA_real_
  lsalpha_init <- NA_real_
  lstheta_init <- NA_real_
  mup_init <- NA_real_
  lsp_init <- NA_real_
  muA_ref_init <- NA_real_
  muB_ref_init <- NA_real_
  
  if(!is.null(init)){
    mup_init     <- init$par.fixed["mu_prev"]
    Lambda_init  <- init$par.fixed["Lambda"]
    Theta_init   <- init$par.fixed["Theta"]
    beta_init    <- init$par.fixed["beta"]
    lsp_init     <- init$par.fixed["log_sigma_prev"]
    lsalpha_init <- init$par.fixed["log_sigma_alpha"]
    lstheta_init <- init$par.fixed["log_sigma_theta"]
    muA_ref_init <- init$par.fixed["mu_A_ref"]
    muB_ref_init <- init$par.fixed["mu_B_ref"]
  }
  
  if (!is.finite(Lambda_init)) {Lambda_init <- 0}
  if (!is.finite(Theta_init)) {Theta_init <- 0}
  if (!is.finite(beta_init)) {beta_init <- 0}
  if (!is.finite(lsalpha_init)) {lsalpha_init <- 0.5*log(0.5)}
  if (!is.finite(lstheta_init)) {lstheta_init <- 0.5*log(0.125)}
  reit <- getREIT(Lambda=Lambda_init,Theta=Theta_init,beta=beta_init,sigma2_alpha=NA,sigma2_theta=NA)
  if (!is.finite(muA_ref_init)) {muA_ref_init <- stats::qlogis(mean(c(stats::plogis(reit$mu_A.sens),0.99)))}
  if (!is.finite(muB_ref_init)) {muB_ref_init <- stats::qlogis(mean(c(stats::plogis(reit$mu_B.spec),0.99)))}
  if (!is.finite(mup_init)) {mup_init <- 0}
  if (!is.finite(lsp_init)) {lsp_init <- log(0.5)}

  # Construct Z and Z_pred
  if(llsub== 1){
    Z <- matrix(1, nrow = nrow(X), ncol = 1)
    Z_pred <- matrix(1, nrow = 1, ncol = 1)
  } else {
    Z   <- model.matrix(~ subgroup, data = X)
    sub <- data.frame(subgroup = factor(lsub, levels = lsub))
    Z_pred <- model.matrix(~ subgroup, data = sub)
  }
  ngroup <- ncol(Z)

  dat2 <- list(
    y11    = X$y11,
    y10    = X$y10,
    y01    = X$y01,
    y00    = X$y00,
    Z      = Z,
    Z_pred = Z_pred,
    spec   = if (is.null(spec)) 0.8 else spec
  )
  
  parameters <- list(
    prev_coef           = c(mup_init,    rep(0,ngroup-1)),
    accuracy_coef       = c(Lambda_init, rep(0,ngroup-1)),
    threshold_coef      = c(Theta_init,  rep(0,ngroup-1)),
    shape_coef          = c(beta_init,   rep(0,ngroup-1)),
    log_sigma_prev_coef = c(lsp_init,rep(0,ngroup-1)), 
    log_sigma_alpha = lsalpha_init,
    log_sigma_theta = lstheta_init,
    mu_A_ref = muA_ref_init,
    mu_B_ref = muB_ref_init,
    
    prevu = rep(0, n_study),
    alpha = rep(0, n_study),
    theta = rep(0, n_study)
  )

  map <- list()
  if("accuracy" %in% constrain){
    map$accuracy_coef = factor(c(1, rep(NA,ngroup-1)))
  }
  if("threshold" %in% constrain){
    map$threshold_coef = factor(c(1, rep(NA,ngroup-1)))
  }
  if("shape" %in% constrain){
    map$shape_coef = factor(c(1, rep(NA,ngroup-1)))
  }
  if("shape_zero" %in% constrain) {
    parameters$shape_coef <- rep(0,ngroup)
    map$shape_coef <- factor(rep(NA,ngroup))
  }
  if("sigma2_alpha" %in% constrain) {
    parameters$log_sigma_alpha <- log(.Machine$double.eps)
    map$log_sigma_alpha <- factor(NA)
  }
  
  if("sigma2_theta" %in% constrain) {
    parameters$log_sigma_theta <- log(.Machine$double.eps)
    map$log_sigma_theta <- factor(NA)
  }
  
  if(prev_variances=="common"){
    map$log_sigma_prev_coef = factor(c(1, rep(NA,ngroup-1)))
  }
  
  dat2$model <- "RutterGatsonisSubgroupLCA"
  ## Objective function ##
  obj <- TMB::MakeADFun(dat2,
                        parameters,
                        map=if(length(map) == 0) NULL else map,
                        random = c("prevu", "alpha", "theta"),
                        silent=!verbose,
                        DLL = "dtametaTMB_TMBExports")
  
  fit  <- stats::nlminb(obj$par, 
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
  
  # Standard errors
  rep  <- TMB::sdreport(obj) 
  rep1 <- summary(rep,select="report")#p.value=TRUE)
  rep2 <- rename_hsrocsubLCA_rows(rep1,lsub)
  
  lspec <- length(dat2$spec)
  ## Get sensitivities and specificities
  qq   <- stats::qnorm(1-(1-conflevel)/2)
  rlse <- which(rownames(rep2)=="logitsens")
  sesp <- data.frame(subgroup=rep(lsub,each=lspec),
                     spec=rep(dat2$spec,llsub),
                     conflevel=conflevel,
                     logitsens=rep2[rlse,"Estimate"],
                     Std_Error=rep2[rlse,"Std. Error"],
                     CI_Lower=NA,
                     CI_Upper=NA)
  sesp$CI_Lower     <- with(sesp,logitsens-qq*Std_Error)
  sesp$CI_Upper     <- with(sesp,logitsens+qq*Std_Error)
  sesp$Sens         <- with(sesp,stats::plogis(logitsens))
  sesp$SensCI_Lower <- with(sesp,stats::plogis(CI_Lower))
  sesp$SensCI_Upper <- with(sesp,stats::plogis(CI_Upper))
  ##
  ### Prevalence
  rlpr           <- paste0("mu_prev.",lsub)
  prre           <- as.data.frame(rep2[c(rlpr,
                                         "mu_A.ref",
                                         "mu_B.ref"),])
  prre$type      <- c(rep("Prev",llsub),"Sens","Spec")
  prre$Orig      <- with(prre,stats::plogis(Estimate))
  prre$conflevel <- conflevel
  prre$CI_Lower  <- with(prre,stats::plogis(Estimate-qq*`Std. Error`))
  prre$CI_Upper  <- with(prre,stats::plogis(Estimate+qq*`Std. Error`))
  prre           <- prre[,c("type","Orig","conflevel","CI_Lower","CI_Upper")]
  colnames(prre) <- c("type","Estimate","conflevel","CI_Lower","CI_Upper")
  ## Recover Reitsma parameters
  lamb <- paste0("Lambda_",lsub)
  thet <- paste0("Theta_",lsub)
  bet  <- paste0("beta_",lsub)
  Lambda <- rep2[lamb,"Estimate"]
  Theta  <- rep2[thet,"Estimate"]
  beta   <- rep2[bet,"Estimate"]
  sigma2_alpha <- rep2["sigma2_alpha","Estimate"]
  sigma2_theta <- rep2["sigma2_theta","Estimate"]
  reit <- reit2 <- c()
  for(i in seq_len(llsub)){
    reit  <- getREIT(Lambda[i], Theta[i], beta[i], sigma2_alpha, sigma2_theta)
    reit2 <- rbind(reit2,reit)
  }
  rownames(reit2) <- lsub

    # Get empirical Bayes estimates
  u_est <- rep$par.random
  u_var <- rep$diag.cov.random
  u     <- data.frame(raneff=names(u_est),
                      u=u_est,
                      u_var=u_var)
  
  names(Lambda) <- sub("^Lambda_", "", names(Lambda))
  names(Theta)  <- sub("^Theta_", "",  names(Theta))
  names(beta)   <- sub("^beta_", "",   names(beta))
  
  theta_i <- u[u$raneff=="theta","u"]
  alpha_i <- u[u$raneff=="alpha","u"]

  var_a   <- u[u$raneff=="alpha","u_var"]
  var_t   <- u[u$raneff=="theta","u_var"]
  
  if ("sigma2_alpha" %in% constrain ||
      "sigma2_theta" %in% constrain) {
    cov_at  <- rep(0, n_study)
  } else {
    rep3    <- TMB::sdreport(obj,getJointPrecision=TRUE)$jointPrecision
    cov_at  <- tryCatch({
      idx_a   <- which(colnames(rep3)=="alpha")
      idx_t   <- which(colnames(rep3)=="theta")
      vcov    <- solve(as.matrix(rep3))[idx_a,idx_t]
      cov_at  <- diag(vcov)
      cov_at}, error = function(e){
        warning(
          "The joint random-effects precision matrix could not be inverted. ",
          "The conditional posterior covariance between the study-specific ",
          "accuracy and threshold effects has been set to zero when ",
          "calculating empirical Bayes sensitivity and specificity variances. ",
          "Study-specific confidence intervals may be approximate. ",
          "Original error: ",
          conditionMessage(e))
        rep(0,n_study)})
  }
  
  
  lsens_study <-  (Theta[X$subgroup] + theta_i + (Lambda[X$subgroup]+alpha_i)/2)*exp(-beta[X$subgroup]/2)
  lspec_study <- -(Theta[X$subgroup] + theta_i - (Lambda[X$subgroup]+alpha_i)/2)*exp( beta[X$subgroup]/2)
  
  X$sens_eb  <- stats::plogis(lsens_study)
  X$spec_eb  <- stats::plogis(lspec_study)
  
  X$lsens_eb_var <- exp(-beta[X$subgroup])*(var_t+1/4*var_a+cov_at)
  X$lspec_eb_var <-  exp(beta[X$subgroup])*(var_t+1/4*var_a-cov_at)
  
  ###
  res <- list(
    data         = X,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = sesp,
    prevref      = prre,
    Reitsma_recovered = reit2,
    constrain    = constrain,
    subgroups    = lsub,
    prev_variances = prev_variances
  )
  class(res) <- c("RutterGatsonisSubgroupLCA","RutterGatsonis","DTAmodel")
  return(res)
}

