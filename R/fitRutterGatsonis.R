#' Fit the Rutter and Gatsonis (HSROC) model
#'
#' Fits the hierarchical summary receiver operating characteristic (HSROC)
#' model as proposed by Rutter and Gatsonis for meta-analysis of diagnostic
#' test accuracy (DTA) studies using Template Model Builder (TMB).
#'
#' @details
#' The function internally transforms the data into long format and fits the model
#' via maximum likelihood using TMB. Random effects are included for study-specific
#' accuracy and threshold parameters.
#'
#' Reitsma parameterization is recovered from the fitted HSROC parameters.
#' 
#' The \code{constrain} argument allows simplified HSROC models to be
#' fitted by fixing selected parameters at zero. These constrained models
#' may improve numerical stability in sparse datasets and can be used to
#' reproduce several simplified models presented in the Cochrane Handbook
#' for Diagnostic Test Accuracy Reviews.
#'
#' Constraints are imposed using parameter mapping within Template Model
#' Builder (TMB), so constrained parameters are excluded from numerical
#' optimization and treated as fixed constants.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param constrain Optional character vector specifying model parameters
#'   that should be fixed at zero during estimation.
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
#'     \item{"shape"}{
#'       Fix the HSROC shape parameter (beta) at zero, resulting in a symmetric
#'       summary ROC curve.
#'     }
#'   }
#'
#'   Multiple constraints may be specified simultaneously, e.g.
#'   \code{constrain = c("sigma2_alpha", "shape")}.
#'
#'   If \code{NULL} (default), the unconstrained HSROC model is fitted.
#' @param spec Optional specificity value at which sensitivity is estimated.
#' If \code{NULL}, the median observed specificity is used as a proxy.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return
#' An object of class \code{"RutterGatsonis"} containing:
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
#' @importFrom stats nlminb median qnorm
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
#' 
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#'
#' @examples
#' data("RF")
#' fit <- fitRutterGatsonis(
#'   data = RF,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study
#' )
#' summary(fit)
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats complete.cases nlminb median qnorm qlogis sd
#' @note Requires a compiled TMB model named \code{"RutterGatsonis"}.
#' @export
fitRutterGatsonis <- function(data,
                              TP,
                              FP,
                              FN,
                              TN,
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
  
  X <- XP <- check_preprocess_data(data,
                                   TP, FP, FN, TN,
                                   study=study,
                                   conflevel=conflevel)
  XP <- getXP(X=XP)
  
  ### Get initial values
  logit_sens   <- stats::qlogis(pmin(pmax(XP$sens,0.005),0.995))
  logit_spec   <- stats::qlogis(pmin(pmax(XP$spec,0.005),0.995))
  muA_init     <- mean(logit_sens)
  muB_init     <- mean(logit_spec)
  sA_init      <- max(stats::sd(logit_sens),1e-05)
  sB_init      <- max(stats::sd(logit_spec),1e-05)
  rAB_init     <- max(min(stats::cor(logit_sens,logit_spec),0.99),-0.99)
  sAB_init     <- rAB_init*sA_init*sB_init
  init <- getRUGA(lspec    = muA_init,
                  lsens    = muB_init,
                  sigma_a  = sA_init,
                  sigma_b  = sB_init,
                  sigma_ab = sAB_init)
  init["sigma2_alpha"] <- max(init["sigma2_alpha"], 1e-10)
  init["sigma2_theta"] <- max(init["sigma2_theta"], 1e-10)
  
  
  ### How do I fit the model?
  n_study <- nrow(X)
  Y <- reshapeX_RUGA(X) 
  
  dat2 <- list(
    y = Y$y,
    n = Y$n,
    x = Y$x,
    spec = if (is.null(spec)) stats::median(XP$spec) else spec,
    study = Y$recordid - 1  # 0-based
  )
  
  parameters <- list(
    Lambda = init$Lambda ,
    Theta = init$Theta ,
    beta = init$beta ,
    log_sigma_alpha = log(sqrt(init$sigma2_alpha)),
    log_sigma_theta = log(sqrt(init$sigma2_theta)),
    alpha = rep(0, n_study),
    theta = rep(0, n_study)
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
  
  dat2$model <- "RutterGatsonis"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("alpha", "theta"),
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
  
  ### Recover Reitsma parameters
  Lambda <- rep$par.fixed["Lambda"]
  Theta  <- rep$par.fixed["Theta"]
  beta   <- rep$par.fixed["beta"]
  sigma2_alpha <- rep$value["sigma2_alpha"]
  sigma2_theta <- rep$value["sigma2_theta"]
  reit   <- getREIT(Lambda, Theta, beta, sigma2_alpha, sigma2_theta)
  # How to get sensitivities
  qq   <- qnorm(1-(1-conflevel)/2)
  rlse <- which(rownames(rep2)=="logitsens")
  sesp <- data.frame(spec=dat2$spec,
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
  # Result object
  res <- list(
    data         = XP,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = sesp,
    Reitsma_recovered = reit,
    constrain    = constrain
  )
  
  # Assign class
  class(res) <- c("RutterGatsonis","Cochrane")
  
  return(res)
}