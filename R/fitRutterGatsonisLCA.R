#' Fit Rutter and Gatsonis LCA Model
#'
#' Fits the Rutter and Gatsonis latent class model for diagnostic test accuracy (DTA),
#' allowing for an imperfect reference standard under conditional independence.
#'
#' @param data A data.frame containing study-level data.
#' @param y11 Number testing positive with both index and reference tests (column name).
#' @param y10 Number testing positive with index test but negative with reference test (column name).
#' @param y01 Number testing negative with index test but positive with reference test (column name).
#' @param y00 Number testing negative with both index and reference tests (column name).
#' @param study Study identifier (column name). 
#' @param constrain Optional character vector specifying model parameters
#'   of the index test that should be fixed at zero during estimation.
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
#'     \item{"shape"}{
#'       Fix the HSROC shape parameter (beta) at zero, resulting in a symmetric
#'       summary ROC curve.
#'     }
#'   }
#'
#'   Multiple constraints may be specified simultaneously, e.g.
#'   \code{constrain = c("sigma2_alpha", "shape")}.
#'
#'   If \code{NULL} (default), the unconstrained HSROC LCA model is fitted.
#' @param spec Optional specificity value at which sensitivity is estimated.
#' If \code{NULL}, a specificity of 0.8 is used as a proxy.
#' 
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' 
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return
#' An object of class \code{"RutterGatsonisLCA"} containing:
#' \describe{
#'   \item{data}{Processed input data with derived quantities.}
#'   \item{fit}{Optimization result from \code{nlminb}.}
#'   \item{sdreport}{TMB standard report.}
#'   \item{sdreport2}{Summary of reported parameters.}
#'   \item{sensspec}{Estimated sensitivity at given specificity with confidence intervals.}
#'   \item{Reitsma_recovered}{Recovered parameters in the Reitsma parameterization.}
#'   \item{constrain}{Parameters fixed at zero.}
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb median qnorm plogis plogis
#'
#' @examples
#' data("anticcp")
#' fit <- fitRutterGatsonisLCA(
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

fitRutterGatsonisLCA <- function(data,
                                 y11,
                                 y10,
                                 y01,
                                 y00,
                                 study,
                                 conflevel=0.95,
                                 constrain=NULL,
                                 spec=NULL,
                                 verbose=FALSE){
  
  allowed_constraints <- c(
    "sigma2_alpha",
    "sigma2_theta",
    "shape"
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
  
  excluded <- !stats::complete.cases(dat)
  if (any(excluded)) {
    removed_studies <- unique(dat$study[excluded])
    message(
      "Removed rows with missing values for studies: ",
      paste(removed_studies, collapse = ", ")
    )
  }
  X <- dat[stats::complete.cases(dat), ]
  X <- check2_data(dat=X,conflevel=conflevel)
  X$n     <- with(X, y11+y10+y01+y00)
  n_study <- nrow(X)
  
  ### Get initial values
  init <- fitRutterGatsonis(data=X,
                            TP=y11,
                            FP=y10,
                            FN=y01,
                            TN=y00,
                            study=study)
  
  Theta        <- init$sdreport2["Theta","Estimate"]
  Lambda       <- init$sdreport2["Lambda","Estimate"]
  beta         <- init$sdreport2["beta","Estimate"]
  sigma2_alpha <- init$sdreport2["sigma2_alpha","Estimate"]
  sigma2_theta <- init$sdreport2["sigma2_theta","Estimate"]
  
  mu_A.ref     <- stats::qlogis(mean(c(stats::plogis(as.numeric(init$Reitsma_recovered["mu_A.sens"])),0.99)))
  mu_B.ref     <- stats::qlogis(mean(c(stats::plogis(as.numeric(init$Reitsma_recovered["mu_B.spec"])),0.99)))
  
  prev_i       <- with(X,(y11+y01+0.5)/(y11+y10+y01+y00+1))
  mu_prev      <- mean(stats::qlogis(prev_i))
  sigma2_prev  <- stats::var(stats::qlogis(prev_i))
  
  prevu        <- rep(0,n_study)
  theta        <- rep(0,n_study)
  alpha        <- rep(0,n_study)
  
  parameters <- list(
    mu_prev = mu_prev,
    Lambda  = Lambda,
    Theta   = Theta,
    beta    = beta,
    
    log_sigma_prev   = 0.5 * log(sigma2_prev),
    log_sigma_alpha  = 0.5 * log(sigma2_alpha),
    log_sigma_theta  = 0.5 * log(sigma2_theta),
    
    mu_A_ref = mu_A.ref,
    mu_B_ref = mu_B.ref,
    
    prevu = prevu,
    theta = theta,
    alpha = alpha
  )
  
  dat2 <- list(
    y11  = X$y11,
    y10  = X$y10,
    y01  = X$y01,
    y00  = X$y00,
    spec = if (is.null(spec)) 0.8 else spec
  )
  
  #################
  ### Constrain ###
  #################
  
  map <- list()
  if ("sigma2_alpha" %in% constrain) {
    parameters$log_sigma_alpha <- log(.Machine$double.eps)
    map$log_sigma_alpha <- factor(NA)
  }
  
  if ("sigma2_theta" %in% constrain) {
    parameters$log_sigma_theta <- log(.Machine$double.eps)
    map$log_sigma_theta <- factor(NA)
  }
  
  if ("shape" %in% constrain) {
    parameters$beta <- 0
    map$beta <- factor(NA)
  }
  
  dat2$model <- "RutterGatsonisLCA"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("prevu","alpha", "theta"),
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
  
  # Recover Reitsma Parameters
  indexREIT <- getREIT(Lambda=rep2["Lambda","Estimate"],
                       Theta=rep2["Theta","Estimate"],
                       beta=rep2["beta","Estimate"],
                       sigma2_alpha = rep2["sigma2_alpha","Estimate"],
                       sigma2_theta = rep2["sigma2_theta","Estimate"])
  
  
  # Get empirical Bayes estimates
  u_est <- rep$par.random
  u_var <- rep$diag.cov.random
  u     <- data.frame(raneff=names(u_est),
                      u=u_est,
                      u_var=u_var)
  
  Theta   <- rep2["Theta","Estimate"] 
  theta_i <- u[u$raneff=="theta","u"]
  Lambda  <- rep2["Lambda","Estimate"]
  alpha_i <- u[u$raneff=="alpha","u"]
  beta    <- rep2["beta","Estimate"]
  
  var_a   <- u[u$raneff=="alpha","u_var"]
  var_t   <- u[u$raneff=="theta","u_var"]
  
  if ("sigma2_alpha" %in% constrain ||
      "sigma2_theta" %in% constrain) {
    cov_at  <- rep(0, n_study)
  } else {
    rep3    <- TMB::sdreport(obj,getJointPrecision=TRUE)$jointPrecision
    idx_a   <- which(colnames(rep3)=="alpha")
    idx_t   <- which(colnames(rep3)=="theta")
    vcov    <- solve(as.matrix(rep3))[idx_a,idx_t]
    cov_at  <- diag(vcov)
  }

  lsens_study <- (Theta + theta_i + (Lambda+alpha_i)/2)*exp(-beta/2)
  lspec_study <- -(Theta + theta_i - (Lambda+alpha_i)/2)*exp(beta/2)
  
  X$sens_eb  <- stats::plogis(lsens_study)
  X$spec_eb  <- stats::plogis(lspec_study)
  
  X$lsens_eb_var <- exp(-beta)*(var_t+1/4*var_a+cov_at)
  X$lspec_eb_var <-  exp(beta)*(var_t+1/4*var_a-cov_at)
  
  # How to get sensitivities
  qq   <- stats::qnorm(1-(1-conflevel)/2)
  ### Sensitivity and Specificity
  rlse <- which(rownames(rep2)=="logitsens")
  sesp <- data.frame(spec=dat2$spec,
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
  # Result object
  res <- list(
    data         = X,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = sesp,
    Reitsma_recovered = indexREIT,
    constrain    = constrain
  )
  
  # Assign class
  class(res) <- c("RutterGatsonisLCA","CochraneLCA")
  
  return(res)
}