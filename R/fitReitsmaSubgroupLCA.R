#' Fit Reitsma Subgroup LCA Model
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
#'   \item \code{sdreport}: fitted model object.
#'   \item \code{sdreport2}: parameter estimates with SE.
#'   \item \code{vcov}: variance-covariance matrix.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{LRDOR}: Diagnostic odds ratio and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{subgroups}: The subgroup levels used in the model fit.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#'   \item \code{variances}: Variance structure used in the fitted model.
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
  
  excluded <- !stats::complete.cases(X)
  if (any(excluded)) {
    removed_studies <- unique(X$study[excluded])
    message(
      "Removed rows with missing values for studies: ",
      paste(removed_studies, collapse = ", ")
    )
  }
  X <- X[stats::complete.cases(X), ]
  
  
  X$subgroup <- factor(X$subgroup)
  lsub       <- levels(X$subgroup)
  n_study    <- nrow(X)
  G          <- length(lsub)
  X$n     <- with(X, y11+y10+y01+y00)

  ### Get initial values
  init <- fitReitsmaLCA(data=X,
                        y11=y11,
                        y10=y10,
                        y01=y01,
                        y00=y00,
                        study=study)

  mu_prev        <- init$sdreport2["mu_prev","Estimate"]
  mu_A.index     <- init$sdreport2["mu_A.index","Estimate"]
  mu_B.index     <- init$sdreport2["mu_B.index","Estimate"]
  sigma2_prev    <- init$sdreport2["sigma2_prev","Estimate"]
  sigma2_A.index <- init$sdreport2["sigma2_A.index","Estimate"]
  sigma2_B.index <- init$sdreport2["sigma2_B.index","Estimate"]
  sigma_AB.index <- init$sdreport2["sigma_AB.index","Estimate"]
  rho_AB.index   <- init$sdreport2["rho_AB.index","Estimate"]
  mu_A.ref       <- init$sdreport2["mu_A.ref","Estimate"]
  mu_B.ref       <- init$sdreport2["mu_B.ref","Estimate"]
  
  parameters <- list(
    mu_prev    = rep(mu_prev,G),
    mu_A_index = rep(mu_A.index,G),
    mu_B_index = rep(mu_B.index,G),
    mu_A_ref   = mu_A.ref,
    mu_B_ref   = mu_B.ref,
    
    log_sigma_prev     = rep(0.5*log(sigma2_prev),G),
    log_sigma_A_index  = rep(0.5*log(sigma2_A.index),G),
    log_sigma_B_index  = rep(0.5*log(sigma2_B.index),G),
    theta_AB_index     = rep(atanh(rho_AB.index),G),
    
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
  
  # Reports
  rep  <- TMB::sdreport(obj)
  rep1 <- summary(rep, select = "report")
  rep2 <- rename_reitsubLCA_rows(rep1,lsub)
  rep2
  
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
  sesp$type[grepl("^mu_prev", rownames(sesp))] <- "Prev"
  sesp$type[grepl("^mu_A",    rownames(sesp))] <- "Sens"
  sesp$type[grepl("^mu_B",    rownames(sesp))] <- "Spec"
  sesp           <- sesp[,c("type","Estimate","conflevel","CI_Lower","CI_Upper")]
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
    sensspec     = sesp,
    LRDOR        = lrdor2,
    RutterGatsonis_recovered = ruga2,
    subgroups    = lsub,
    constrain    = constrain,
    variances    = variances
  )
  
  # Assign class
  class(res) <- c("ReitsmaSubgroupLCA")
  
  return(res)
}