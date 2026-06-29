#' Fit the Rutter and Gatsonis (HSROC) model
#'
#' Fits the hierarchical summary receiver operating characteristic (HSROC)
#' model as proposed by Rutter and Gatsonis for meta-analysis of diagnostic
#' test accuracy (DTA) studies using Template Model Builder (TMB).
#'
#' @description
#' This function estimates the HSROC model parameters based on study-level
#' 2x2 table data (true positives, false positives, false negatives, true negatives).
#'
#' @details
#' The function internally transforms the data into long format and fits the model
#' via maximum likelihood using TMB. Random effects are included for study-specific
#' accuracy and threshold parameters.
#'
#' Reitsma parameterization is recovered from the fitted HSROC parameters.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param spec Optional specificity value at which sensitivity is estimated.
#' If \code{NA}, the median observed sensitivity is used as a proxy.
#' @param verbose Logical. Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return
#' An object of class \code{"RutterGatsonis"} containing:
#' \describe{
#'   \item{data}{Processed input data with derived quantities.}
#'   \item{fit}{Optimization result from \code{nlminb}.}
#'   \item{sdreport}{TMB standard report.}
#'   \item{sdreport2}{Summary of reported parameters.}
#'   \item{specsens}{Estimated sensitivity at given specificity with confidence intervals.}
#'   \item{Reitsma_recovered}{Recovered parameters in the Reitsma parameterization.}
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
#' @seealso
#' \code{\link{fitRutterGatsonis}} (print/summary methods)
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
fitRutterGatsonis <- function(data,TP,FP,FN,TN,study,conflevel=0.95,spec=NA,verbose=FALSE){
  
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
      FN = data[[FN_col]])
    
  excluded <- !stats::complete.cases(dat)
  if (any(excluded)) {
    removed_studies <- unique(dat$study[excluded])
    message(
      "Removed rows with missing values for studies: ",
      paste(removed_studies, collapse = ", ")
    )
  }
  
  dat <- dat[stats::complete.cases(dat), ]
  # Validation
  numeric_cols <- c("TP", "TN", "FP", "FN")
  non_numeric <- numeric_cols[!sapply(dat[numeric_cols], is.numeric)]
  if (length(non_numeric) > 0) {
    stop("Columns must be numeric: ", paste(non_numeric, collapse = ", "))
  }
  
  # Columns
  count_cols <- c("TP", "TN", "FP", "FN")
  # Check for non-integers or negative values
  invalid_counts <- sapply(dat[count_cols], function(x) {
    any(x < 0 | x != floor(x), na.rm = TRUE)
  })
  
  if (any(invalid_counts)) {
    stop("Columns TP, TN, FP, FN must contain non-negative integer counts.")
  }
  
  X <- XP <- dat
  XP$n1    <- XP$TP+XP$FN
  XP$n0    <- XP$FP+XP$TN
  XP$sens  <- XP$TP / XP$n1
  XP$spec  <- XP$TN / XP$n0
  XP$recordid <- 1:nrow(X)
  
  ### Get initial values
  logit_sens   <- stats::qlogis(pmin(pmax(XP$sens,0.005),0.995))
  logit_spec   <- stats::qlogis(pmin(pmax(XP$spec,0.005),0.995))
  muA_init     <- mean(logit_sens)
  muB_init     <- mean(logit_spec)
  sA_init      <- max(stats::sd(logit_sens),1e-05)
  sB_init      <- max(stats::sd(logit_spec),1e-05)
  rAB_init     <- max(min(stats::cor(logit_sens,logit_spec),0.99),-0.99)
  sAB_init     <- rAB_init*sA_init*sB_init
  Lambda_init       <- (((sB_init/sA_init)**0.5) * muA_init) + ((sA_init/sB_init)**0.5 *muB_init)
  Theta_init        <- 0.5*((((sB_init/sA_init)**0.5 )*muA_init) - (((sA_init/sB_init)**0.5) *muB_init))
  beta_init         <- log(sB_init/sA_init)
  sigma2_alpha_init <- 2*((sA_init*sB_init) + sAB_init)
  sigma2_alpha_init <- max(sigma2_alpha_init, 1e-10)
  sigma2_theta_init <- 0.5*((sA_init*sB_init) - sAB_init)
  sigma2_theta_init <- max(sigma2_theta_init, 1e-10)
  
  
  ### How do I fit the model?
  X$y1 <- X$TP
  X$y0 <- X$FP  
  X$n1 <- X$TP+X$FN
  X$n0 <- X$FP+X$TN
  n_study <- nrow(X)
  
  ### Reshape the data from wide to long format. ###
  Y = reshape(X, direction="long", varying=list(c("n1", "n0"), c("y1","y0")),
              timevar="x", times=c(0.5,-0.5), v.names=c("n","y")) 
  Y = Y[order(Y$id),]  
  
  dat2 <- list(
    y = Y$y,
    n = Y$n,
    x = Y$x,
    spec = ifelse(is.na(spec),
                  stats::median(XP$spec),
                  spec),
    study = Y$id - 1  # 0-based
  )
  
  parameters <- list(
    Lambda = Lambda_init,
    Theta = Theta_init,
    beta = beta_init,
    log_sigma_alpha = log(sqrt(sigma2_alpha_init)),
    log_sigma_theta = log(sqrt(sigma2_theta_init)),
    alpha = rep(0, n_study),
    theta = rep(0, n_study)
  )
  
  dat2$model <- "RutterGatsonis"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
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
  b      <- exp(beta/2)
  sigma2_alpha <- rep$value["sigma2_alpha"]
  sigma2_theta <- rep$value["sigma2_theta"]
  
  rep5 <- data.frame(
    mu_A.sens = b**-1*(Theta+0.5*Lambda),
    mu_B.spec = -b*(Theta-0.5*Lambda),
    sigma2_A.sens = b**-2*(sigma2_theta+0.25*sigma2_alpha),
    sigma2_B.spec = b**2*(sigma2_theta+0.25*sigma2_alpha),
    sigma_AB = -(sigma2_theta-0.25*sigma2_alpha),
    row.names="Estimate (recovered)")
  # How to get sensitivities
  qq   <- qnorm(1-(1-conflevel)/2)
  rlse <- which(rownames(rep2)=="logitsens")
  rep6 <- data.frame(spec=dat2$spec,
                     conflevel=conflevel,
                     logitsens=rep2[rlse,"Estimate"],
                     Std_Error=rep2[rlse,"Std. Error"],
                     CI_Lower=NA,
                     CI_Upper=NA)
  rep6$CI_Lower   <- with(rep6,logitsens-qq*Std_Error)
  rep6$CI_Upper   <- with(rep6,logitsens+qq*Std_Error)
  rep6$Sens       <- with(rep6,plogis(logitsens))
  rep6$SensCI_Lower <- with(rep6,plogis(CI_Lower))
  rep6$SensCI_Upper <- with(rep6,plogis(CI_Upper))
  rep6
  # Result object
  res <- list(
    data         = XP,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = rep6,
    Reitsma_recovered = rep5
  )
  
  # Assign class
  class(res) <- c("RutterGatsonis","Cochrane")
  
  return(res)
}