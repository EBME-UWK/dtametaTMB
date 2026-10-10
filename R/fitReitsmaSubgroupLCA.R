#' Fit the Reitsma Subgroup LCA Model
#'
#' Fits the Reitsma latent class model with a single categorical covariate,
#' allowing for an imperfect reference standard under conditional independence.
#'
#' @param data A data.frame containing study-level data.
#' @param y11 Number testing positive with both index and reference tests (column name).
#' @param y10 Number testing positive with index test but negative with reference test (column name).
#' @param y01 Number testing negative with index test but positive with reference test (column name).
#' @param y00 Number testing negative with both index and reference tests (column name).
#' @param study Study identifier (column name).
#' @param subgroup A single categorical study-level subgroup variable (column name).
#' @param constrain Optional character string specifying a simplified
#'   covariance structure for the index test in the Reitsma LCA model.
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
#'   These random-effects simplifications are currently only available for 
#'   \code{variances="common"}.
#'   
#' @param sensspec_constrain Optional character vector specifying
#'   restrictions on subgroup-specific logit-sensitivity and/or
#'   logit-specificity parameters for the index test.
#'
#'   This can be useful for testing whether subgroup differences
#'   are present in sensitivity, specificity, or both.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{"sens"}}{
#'       Constrain all subgroup-specific logit-sensitivity parameters
#'       to be equal. Subgroup differences are therefore only allowed
#'       in logit-specificity.
#'     }
#'     \item{\code{"spec"}}{
#'       Constrain all subgroup-specific logit-specificity parameters
#'       to be equal. Subgroup differences are therefore only allowed
#'       in logit-sensitivity.
#'     }
#'   }
#'
#'   Both constraints may be specified simultaneously, e.g.
#'   \code{sensspec_constrain = c("sens", "spec")}, which forces
#'   all subgroup-specific sensitivity and specificity parameters
#'   to be equal across subgroups.
#'
#'   If \code{NULL} (default), separate sensitivity and specificity
#'   parameters are estimated for each subgroup.
#'   
#' @param variances Whether the between-study random-effects variance-covariance 
#'   matrix of the index test should be assumed to be \code{"common"} (default) or \code{"unequal"} 
#'   across subgroups. If \code{"common"}, a single between-study variance-covariance 
#'   matrix is estimated and shared across all subgroups. If \code{"unequal"}, 
#'   subgroup-specific variance-covariance matrices are estimated.
#'
#' @param prev_variances Whether the between-study random-effects variance for 
#'   logit prevalence should be assumed to be \code{"common"} (default) or \code{"unequal"} 
#'   across subgroups. If \code{"common"}, a single between-study variance  is estimated and shared across all subgroups. 
#'   If \code{"unequal"}, subgroup-specific variances are estimated.
#'   
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return A list of class \code{"ReitsmaSubgroupLCA"} with components:
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
#'   \item \code{subgroups}: The subgroup levels used in the model fit.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#'   \item \code{sensspec_constrain}: Whether sensitivities/specificities are assumed common across subgroups.
#'   \item \code{variances}: Variance structure used in the fitted model.
#'   \item \code{prev_variances}: Prevalence variance structure used in the fitted model.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsmaSubgroupLCA(
#'   data = anticcp,
#'   y11 = TP,
#'   y10 = FP,
#'   y01 = FN,
#'   y00 = TN,
#'   study = study,
#'   subgroup = generation
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
#' @importFrom stats nlminb complete.cases qnorm plogis qlogis
#' @export
fitReitsmaSubgroupLCA <- function(data,
                                  y11, y10, y01, y00,
                                  study,
                                  subgroup,
                                  constrain=NULL,
                                  sensspec_constrain=NULL,
                                  variances=c("common","unequal"),
                                  prev_variances=c("common","unequal"),
                                  conflevel=0.95,
                                  verbose=FALSE) {
  
  variances      <- match.arg(variances)
  prev_variances <- match.arg(prev_variances)
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  
  allowed_constraints <- c(
    "sigma_AB.index",
    "sigma2_A.index",
    "sigma2_B.index",
    "all"
  )
  
  if (!is.null(constrain)) {
    if (!is.character(constrain) ||
        length(constrain) != 1L ||
        !constrain %in% allowed_constraints) {
      stop(
        "'constrain' must be one of: ",
        paste(shQuote(allowed_constraints), collapse = ", "),
        " or NULL."
      )
    }
  }
  
  if(variances=="unequal" &&
     !is.null(constrain)) {
    stop(
      "'constrain' is currently not supported when ",
      "variances='unequal'."
    )
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
  G          <- length(lsub)
  X$n        <- with(X, y11+y10+y01+y00)

  ### Get initial values
  init <- tryCatch(fitReitsmaLCA(data=X,
                                 y11=y11,
                                 y10=y10,
                                 y01=y01,
                                 y00=y00,
                                 study=study,
                                 conflevel=conflevel,
                                 constrain=NULL)$sdreport,
                   error = function(e) NULL)

  muA_init <- NA_real_
  muB_init <- NA_real_
  lsA_init <- NA_real_
  lsB_init <- NA_real_
  thetaAB_init <- NA_real_
  mup_init <- NA_real_
  lsp_init <- NA_real_
  muA_ref_init <- NA_real_
  muB_ref_init <- NA_real_
  
  if(!is.null(init)){
    mup_init      <- init$par.fixed["mu_prev"]
    muA_init      <- init$par.fixed["mu_A_index"]
    muB_init      <- init$par.fixed["mu_B_index"]
    lsp_init      <- init$par.fixed["log_sigma_prev"]
    lsA_init      <- init$par.fixed["log_sigma_A_index"]
    lsB_init      <- init$par.fixed["log_sigma_B_index"]
    thetaAB_init  <- init$par.fixed["theta_AB_index"]
    muA_ref_init  <- init$par.fixed["mu_A_ref"]
    muB_ref_init  <- init$par.fixed["mu_B_ref"]
  }
  
  if (!is.finite(muA_init)) {muA_init <- 0}
  if (!is.finite(muB_init)) {muB_init <- 0}
  if (!is.finite(lsA_init)) {lsA_init <- log(0.5)}
  if (!is.finite(lsB_init)) {lsB_init <- log(0.5)}
  if (!is.finite(thetaAB_init)) {thetaAB_init <- 0}
  if (!is.finite(muA_ref_init)) {muA_ref_init <- stats::qlogis(mean(c(stats::plogis(muA_init),0.99)))}
  if (!is.finite(muB_ref_init)) {muB_ref_init <- stats::qlogis(mean(c(stats::plogis(muB_init),0.99)))}
  if (!is.finite(mup_init)) {mup_init <- 0}
  if (!is.finite(lsp_init)) {lsp_init <- log(0.5)}

  parameters <- list(
    mu_prev    = rep(mup_init,G),
    mu_A_index = rep(muA_init,G),
    mu_B_index = rep(muB_init,G),
    mu_A_ref   = muA_ref_init,
    mu_B_ref   = muB_ref_init,
    
    log_sigma_prev     = rep(lsp_init,G),
    log_sigma_A_index  = rep(lsA_init,G),
    log_sigma_B_index  = rep(lsB_init,G),
    theta_AB_index     = rep(thetaAB_init,G),
    
    prevu = rep(0,n_study),
    sensu = rep(0,n_study),
    specu = rep(0,n_study)
  )
  
  dat2 <- list(
    y11 = X$y11,
    y10 = X$y10,
    y01 = X$y01,
    y00 = X$y00,
    group = as.numeric(X$subgroup)-1
  )
  
  ## Constraints ###
  map <- list()
  
  if(variances=="common"){
    map$log_sigma_A_index <- factor(rep(1,G))
    map$log_sigma_B_index <- factor(rep(1,G))
    map$theta_AB_index    <- factor(rep(1,G))
  }
  
  if(variances=="common"){
    if(!is.null(constrain)){
      if (constrain == "sigma_AB.index") {
        parameters$theta_AB_index    <- rep(0,G)
        map$theta_AB_index           <- rep(factor(NA),G)
      }
      if (constrain == "sigma2_A.index") {
        parameters$log_sigma_A_index <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB_index    <- rep(0,G)
        map$log_sigma_A_index        <- rep(factor(NA),G)
        map$theta_AB_index           <- rep(factor(NA),G)
      }
      if (constrain == "sigma2_B.index") {
        parameters$log_sigma_B_index <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB_index    <- rep(0,G)
        map$log_sigma_B_index        <- rep(factor(NA),G)
        map$theta_AB_index           <- rep(factor(NA),G)
      }
      if (constrain == "all") {
        parameters$log_sigma_A_index <- rep(log(.Machine$double.eps),G)
        parameters$log_sigma_B_index <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB_index    <- rep(0,G)
        map$log_sigma_A_index        <- rep(factor(NA),G)
        map$log_sigma_B_index        <- rep(factor(NA),G)
        map$theta_AB_index           <- rep(factor(NA),G)
      }
    }
  }
  
  allowed_sensspec_constraints <- c(
    "sens",
    "spec"
  )
  
  if (!is.null(sensspec_constrain)) {
    if (!is.character(sensspec_constrain)) {
      stop("'sensspec_constrain' must be a character vector or NULL.")
    }
    invalid_constraints <- setdiff(
      sensspec_constrain,
      allowed_sensspec_constraints
    )
    if(length(invalid_constraints) > 0){
      stop(
        "Unknown sensspec constraint(s): ",
        paste(invalid_constraints, collapse=", ")
      )
    }
  }
  
  if (!is.null(sensspec_constrain)) {
    if ("sens" %in% sensspec_constrain) {
      map$mu_A_index <- factor(rep(1,G))
    }
    if ("spec" %in% sensspec_constrain) {
      map$mu_B_index <- factor(rep(1,G))
    }
  }
  
  if(prev_variances=="common"){
      map$log_sigma_prev    <- factor(rep(1,G))
  }
  
  dat2$model = "ReitsmaSubgroupLCA"
  
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
  if (!sdreport(obj)$pdHess) {
    warning(
      "The Hessian matrix is not positive definite. ",
      "Parameter estimates and standard errors may be unreliable. ",
      "Consider checking model specification, data quality, or potential identifiability issues."
    )
  }
  
  # Reports
  rep  <- TMB::sdreport(obj)
  rep1 <- summary(rep, select = "report")
  rep2 <- rename_reitsubLCA_rows(rep1,lsub)

  # Variance covariance matrix of fixed effects
  vcov <- rep$cov
  colnames(vcov) <- rownames(vcov) <- rownames(rep2)
  
  
  # Get empirical Bayes estimates
  u_est <- rep$par.random
  u_var <- rep$diag.cov.random
  u     <- data.frame(raneff=names(u_est),
                      u=u_est,
                      u_var=u_var)
  
  mu_A        <- rep2[grep("^mu_A\\.index\\.", rownames(rep2)), "Estimate"]
  names(mu_A) <- sub("^mu_A\\.index\\.", "", names(mu_A))
  mu_B        <- rep2[grep("^mu_B\\.index\\.", rownames(rep2)), "Estimate"]
  names(mu_B) <- sub("^mu_B\\.index\\.", "", names(mu_B))
  
  lsens_study <- mu_A[X$subgroup] + u[u$raneff=="sensu","u"]
  lspec_study <- mu_B[X$subgroup] + u[u$raneff=="specu","u"]
  
  X$sens_eb  <- stats::plogis(lsens_study)
  X$spec_eb  <- stats::plogis(lspec_study)
  
  X$lsens_eb_var <- u[u$raneff=="sensu","u_var"]
  X$lspec_eb_var <- u[u$raneff=="specu","u_var"]
  
  ### Sensitivity and Specificity
  qq             <- stats::qnorm(1-(1-conflevel)/2)
  sesp           <- as.data.frame(rep2[grepl("^mu_",rownames(rep2)),,drop=FALSE])
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,c("Orig","conflevel","CI_Lower","CI_Upper")]
  colnames(sesp) <- c("Estimate","conflevel","CI_Lower","CI_Upper")
  sesp$type <- NA_character_
  rn             <- rownames(sesp)
  sesp$type[grepl("^mu_prev", rn)] <- "Prev"
  sesp$type[grepl("^mu_A",    rn)] <- "Sens"
  sesp$type[grepl("^mu_B",    rn)] <- "Spec"
  sesp           <- sesp[,c("type","Estimate","conflevel","CI_Lower","CI_Upper")]
  sesp1          <- sesp[grepl("^mu_[AB]\\.index",rn),]
  prevref        <- sesp[grepl("^mu_prev|^mu_A\\.ref|^mu_B\\.ref",rn),]
  ### Diagnostic odds ratios and Likelihood ratios
  lrdor2 <- data.frame()
  for(i in seq_along(lsub)){
    sg      <- lsub[i]
    mu_A.sg <- paste0("mu_A.index.",sg)
    mu_B.sg <- paste0("mu_B.index.",sg)
    lsens   <- rep2[mu_A.sg,"Estimate"]
    lspec   <- rep2[mu_B.sg,"Estimate"]
    S       <- vcov[c(mu_A.sg,mu_B.sg),c(mu_A.sg,mu_B.sg)]
    lrdor   <- getLRDOR(lsens=lsens, lspec=lspec, S=S, conflevel=conflevel)
    rownames(lrdor) <- paste0(lsub[i],": ",rownames(lrdor))
    lrdor2  <- rbind(lrdor2,lrdor)
  }
  ### Recover HSROC parameters
  ruga2 <- data.frame()
  for(i in seq_along(lsub)){
    sg      <- lsub[i]
    mu_A.sg <- paste0("mu_A.index.",sg)
    mu_B.sg <- paste0("mu_B.index.",sg)
    s2_A.sg <- paste0("sigma2_A.index.",sg)
    s2_B.sg <- paste0("sigma2_B.index.",sg)
    s_AB.sg <- paste0("sigma_AB.index.",sg)
    ruga <- getRUGA(lsens=rep2[mu_A.sg,"Estimate"],
                    lspec=rep2[mu_B.sg,"Estimate"],
                    sigma_a=sqrt(rep2[s2_A.sg,"Estimate"]),
                    sigma_b=sqrt(rep2[s2_B.sg,"Estimate"]),
                    sigma_ab=rep2[s_AB.sg,"Estimate"])
    ruga2 <- rbind(ruga2,ruga)
  }
  rownames(ruga2) <- lsub
  
  ##
  res <- list(
    data         = X,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    vcov         = vcov,
    sensspec     = sesp1,
    prevref      = prevref,
    LRDOR        = lrdor2,
    RutterGatsonis_recovered = ruga2,
    subgroups    = lsub,
    constrain    = constrain,
    sensspec_constrain = sensspec_constrain,
    variances    = variances,
    prev_variances = prev_variances
  )
  
  # Assign class
  class(res) <- c("ReitsmaSubgroupLCA","Reitsma","DTAmodel")
  
  return(res)
}