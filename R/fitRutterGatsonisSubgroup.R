#' Fit the Rutter and Gatsonis Subgroup Model
#'
#' The function fits an HSROC model with a single categorical study-level
#' covariate (subgroup). 
#'
#' @details
#' The function fits an HSROC model with a single categorical study-level
#' covariate (subgroup). Separate subgroup effects are estimated for the
#' accuracy and threshold parameters. Optionally, subgroup-specific effects
#' can also be estimated for the shape parameter.
#'
#' The fitted model returns subgroup-specific summary sensitivity estimates
#' evaluated at user-specified specificity values. If \code{spec = NA}, the
#' median observed specificity is used as a proxy.
#'
#' This function is intended as a convenient wrapper for subgroup analyses
#' with a single categorical covariate and provides output suitable for
#' dedicated \code{summary()}, \code{plot()}, and \code{forest()} methods.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param subgroup A single categorical study-level subgroup variable (column name).
#' @param constrain Optional character vector specifying model parameters
#'   that should be constrained during estimation.
#'
#'   This can be useful for sparse data, small meta-analyses, or for
#'   reproducing simplified HSROC models described in the Cochrane
#'   Handbook for Diagnostic Test Accuracy Reviews.
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
#' @param spec Optional specificity value or vector of specificity values at which
#' sensitivity is estimated. If \code{NULL}, the median observed specificity is
#' used as a proxy.
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return
#' An object of class \code{"RutterGatsonisSubgroup"} containing:
#' \describe{
#'   \item{data}{Processed input data with derived quantities.}
#'   \item{fit}{Optimization result from \code{nlminb}.}
#'   \item{sdreport}{TMB standard report.}
#'   \item{sdreport2}{Summary of reported subgroup-specific parameters.}
#'   \item{sensspec}{Estimated subgroup-specific sensitivities at the given
#'   specificity value(s), with confidence intervals.}
#'   \item{subgroups}{The subgroup levels used in the model fit.}
#'   \item{constrain}{Constraints on parameters applied during model fitting.}
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats complete.cases median nlminb qnorm
#'
#' @references
#' Reitsma, J. B., et al. (2005).
#' Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews.
#' \emph{Journal of Clinical Epidemiology}, 58(10), 982--990.
#' \doi{10.1016/j.jclinepi.2005.02.022}
#'
#' Rutter, C. M., & Gatsonis, C. A. (2001).
#' A hierarchical regression approach to meta-analysis of diagnostic test accuracy evaluations.
#' \emph{Statistics in Medicine}, 20(19), 2865--2884.
#' \doi{10.1002/sim.942}
#'
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#' @examples
#' data("RF")
#' fit <- fitRutterGatsonisSubgroup(
#'   data = RF,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study,
#'   subgroup = method
#' )
#' summary(fit)
#'
#' @note Requires a compiled TMB model named \code{"RutterGatsonisReg"}.
#' @export
fitRutterGatsonisSubgroup <- function(data,
                                      TP, FP, FN, TN,
                                      study,
                                      subgroup,
                                      constrain=NULL,
                                      spec=NULL,
                                      conflevel=0.95,
                                      verbose=FALSE){
  
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
  
  TP_col       <- deparse(substitute(TP))
  FP_col       <- deparse(substitute(FP))
  FN_col       <- deparse(substitute(FN))
  TN_col       <- deparse(substitute(TN))
  study_col    <- deparse(substitute(study))
  subgroup_col <- deparse(substitute(subgroup))
  
  X <- data.frame(study = data[[study_col]],
                  TP = data[[TP_col]],
                  TN = data[[TN_col]],
                  FP = data[[FP_col]],
                  FN = data[[FN_col]],
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
  
  XP <- getXP(X=X)
  # Get starting values
  init <- fitRutterGatsonis(data=X,
                            TP=TP,
                            TN=TN,
                            FN=FN,
                            FP=FP,
                            study=study,
                            constrain=NULL)$sdreport
  Lambda_init  <- init$par.fixed["Lambda"]
  Theta_init   <- init$par.fixed["Theta"]
  beta_init    <- init$par.fixed["beta"]
  lsalpha_init <- init$par.fixed["log_sigma_alpha"]
  lstheta_init <- init$par.fixed["log_sigma_theta"]
  # Reshape data
  n_study <- nrow(X)
  Y       <- reshapeX_RUGA(X)
  # Construct Z and Z_pred
  lsub   <- levels(Y$subgroup)
  llsub  <- length(lsub)
  if(llsub== 1){
    Z <- matrix(1, nrow = nrow(Y), ncol = 1)
    Z_pred <- matrix(1, nrow = 1, ncol = 1)
  } else {
    Z   <- model.matrix(~ subgroup, data = Y)
    sub <- data.frame(subgroup = factor(lsub, levels = lsub))
    Z_pred <- model.matrix(~ subgroup, data = sub)
  }
  ngroup <- ncol(Z)
  
  dat2 <- list(
    y      = Y$y,
    n      = Y$n,
    x      = Y$x,
    Z      = Z,
    Z_pred = Z_pred,
    study  = Y$recordid - 1,  # 0-based
    spec   = if (is.null(spec)) stats::median(XP$spec) else spec
  )
  parameters <- list(
    accuracy_coef   = c(Lambda_init,rep(0,ngroup-1)),
    threshold_coef  = c(Theta_init, rep(0,ngroup-1)),
    shape_coef      = c(beta_init,  rep(0,ngroup-1)),
    log_sigma_alpha = lsalpha_init,
    log_sigma_theta = lstheta_init,
    alpha = rep(0, n_study),
    theta = rep(0, n_study)
  )
  
  dat2$model <- "RutterGatsonisReg"
  
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
 
  obj <- TMB::MakeADFun(dat2,
                          parameters,
                          map=if(length(map) == 0) NULL else map,
                          random = c("alpha", "theta"),
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
  rep2 <- rename_hsroc_rows(rep1,lsub)
  
  lspec <- length(dat2$spec)
  ## Get sensitivities and specificities
  qq   <- qnorm(1-(1-conflevel)/2)
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
  sesp$Sens         <- with(sesp,plogis(logitsens))
  sesp$SensCI_Lower <- with(sesp,plogis(CI_Lower))
  sesp$SensCI_Upper <- with(sesp,plogis(CI_Upper))
  sesp
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
  for(i in 1:llsub){
    reit  <- getREIT(Lambda[i], Theta[i], beta[i], sigma2_alpha, sigma2_theta)
    reit2 <- rbind(reit2,reit)
  }
  rownames(reit2) <- lsub
  ###
  res <- list(
    data         = XP,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = sesp,
    Reitsma_recovered = reit2,
    constrain    = constrain,
    subgroups    = lsub
  )
  class(res) <- c("RutterGatsonisSubgroup","CochraneSubgroup")
  return(res)
}

