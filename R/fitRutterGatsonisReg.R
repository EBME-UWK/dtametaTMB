#' Fit the Rutter and Gatsonis Regression Model
#'
#' Fits a general HSROC meta-regression model using user-supplied design
#' matrices for estimation and prediction.
#'
#' This function is intended for advanced users who wish to specify design
#' matrices directly. The user is responsible for constructing appropriate
#' design matrices, choosing prediction covariate patterns through
#' \code{Z_pred}, and interpreting or visualizing the resulting regression
#' model.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param init Optional list of initial parameter values.
#' If \code{NULL} (default), all fixed and random effects are initialized
#' to zero. Advanced users may provide a named list containing:
#' \describe{
#'   \item{accuracy_coef}{Initial values for the accuracy regression coefficients.}
#'   \item{threshold_coef}{Initial values for the threshold regression coefficients.}
#'   \item{shape_coef}{Initial values for the shape regression coefficients.}
#'   \item{log_sigma_alpha}{Initial value for the log standard deviation of the
#'   accuracy random effects.}
#'   \item{log_sigma_theta}{Initial value for the log standard deviation of the
#'   threshold random effects.}
#'   \item{alpha}{Initial values for the study-specific accuracy random effects.}
#'   \item{theta}{Initial values for the study-specific threshold random effects.}
#' }
#' The lengths of \code{accuracy_coef}, \code{threshold_coef}, and
#' \code{shape_coef} must equal the number of columns in \code{Z}.
#' The lengths of \code{alpha} and \code{theta} must equal the number
#' of studies.
#'
#' @param Z Design matrix containing covariate values for model fitting.
#' The current implementation is intended primarily for study-level
#' covariates. In the long-format representation, this typically means
#' that each study is represented by two consecutive rows with identical
#' covariate values. The number of rows of \code{Z} must equal twice the
#' number of studies.
#'
#' @param Z_pred Prediction design matrix used to obtain covariate-specific
#' parameter estimates and SROC curves. Each row defines a covariate pattern
#' at which predictions are evaluated. The number of columns of
#' \code{Z_pred} must equal the number of columns of \code{Z}. If
#' \code{NULL} (default), a single row of zeros is used. For meaningful
#' predictions, users will typically want to specify \code{Z_pred}
#' explicitly.
#'
#' @param map Optional named list of parameter mappings passed to
#'   \code{\link[TMB]{MakeADFun}}.
#'
#'   Parameter mapping allows selected model parameters to be fixed or
#'   constrained during estimation. This can be useful for fitting
#'   simplified HSROC models, imposing equality constraints, or
#'   reproducing model specifications described in the literature.
#'
#'   Each component of \code{map} should be a factor vector with the same
#'   length as the corresponding parameter. Parameters assigned
#'   \code{NA} levels are fixed at their initial values supplied via the
#'   \code{init} argument, whereas parameters sharing the same
#'   factor level are estimated as equal.
#'
#'   This is an advanced feature intended primarily for users familiar
#'   with Template Model Builder (TMB).
#'
#'   If \code{NULL} (default), all model parameters are estimated freely.
#'   
#'   Example \code{map = list(shape_coef=factor(c(1, rep(NA, ncol(Z) - 1))))}.
#'
#' @param spec Optional specificity value or vector of specificity values
#' at which sensitivity is evaluated for each covariate pattern specified
#' in \code{Z_pred}. If \code{NULL}, the median observed specificity is used.
#'
#' @param conflevel Confidence level for confidence intervals.
#' Default is \code{0.95}.
#'
#' @param verbose Logical indicating whether TMB optimization output should
#' be printed (default: \code{FALSE}).
#'
#' @return
#' An object of class \code{"RutterGatsonisReg"} containing:
#' \describe{
#'   \item{data}{Processed input data with derived quantities.}
#'   \item{fit}{Optimization result returned by \code{nlminb}.}
#'   \item{sdreport}{TMB standard report object.}
#'   \item{sdreport2}{Summary of reported model parameters and derived quantities.}
#'   \item{sensspec}{Estimated sensitivities at the specified specificity
#'   value(s), including confidence intervals.}
#'   \item{constrain}{Constraints on parameters applied during model fitting.}
#' }
#'
#' @details
#' By default, all fixed and random effects are initialized at zero.
#' For difficult optimization problems, user-supplied starting values may be
#' provided through the \code{init} argument.
#'
#' The fitted model estimates covariate effects on the HSROC accuracy,
#' threshold, and optionally shape parameters. Predicted sensitivities are
#' obtained by evaluating the resulting SROC curve(s) at the specificity
#' value(s) supplied through \code{spec} and the covariate patterns defined
#' by \code{Z_pred}.
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats complete.cases median nlminb qnorm plogis
#'
#' @references
#' Reitsma, J. B., Glas, A. S., Rutjes, A. W. S., Scholten, R. J. P. M.,
#' Bossuyt, P. M., & Zwinderman, A. H. (2005).
#' Bivariate analysis of sensitivity and specificity produces informative
#' summary measures in diagnostic reviews.
#' \emph{Journal of Clinical Epidemiology}, 58(10), 982--990.
#' \doi{10.1016/j.jclinepi.2005.02.022}
#'
#' Rutter, C. M., & Gatsonis, C. A. (2001).
#' A hierarchical regression approach to meta-analysis of diagnostic test
#' accuracy evaluations.
#' \emph{Statistics in Medicine}, 20(19), 2865--2884.
#' \doi{10.1002/sim.942}
#'
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., &
#' Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#' @examples
#' data("RF")
#'
#' Z <- model.matrix(~ method, data = RF)
#' Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
#'
#' fit <- fitRutterGatsonisReg(
#'   data = RF,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study,
#'   Z = Z
#' )
#'
#' summary(fit)
#'
#' @note Requires a compiled TMB model named \code{"RutterGatsonisReg"}.
#'
#' @export
fitRutterGatsonisReg <- function(data,
                                 TP, FP, FN, TN,
                                 study,
                                 Z,
                                 Z_pred=NULL,
                                 map=NULL,
                                 init=NULL,
                                 spec=NULL,
                                 conflevel=0.95,
                                 verbose=FALSE){
  
  # Construct Z and Z_pred
  ngroup <- ncol(Z)
  if(is.null(Z_pred)){Z_pred <- matrix(0,ncol=ngroup,nrow=1)}
  
  if (ncol(Z) != ncol(Z_pred)) {
    stop("'Z' and 'Z_pred' must have the same number of columns.")
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
    FN = data[[FN_col]])
  
  X <- XP <- check_data(dat,
                        conflevel=conflevel)
  
  XP <- getXP(X)
  n_study <- nrow(X)
  ###
  Y     <- reshapeX_RUGA(X)

  if (nrow(Z) != nrow(Y)) {
    stop(
      "'Z' must contain exactly two rows per study",
      "(one diseased and one non-diseased row)."
    )
  }
  
  dat2 <- list(
    y      = Y$y,
    n      = Y$n,
    x      = Y$x,
    Z      = Z,
    Z_pred = Z_pred,
    study  = Y$recordid - 1,  # 0-based
    spec   = if (is.null(spec)) stats::median(XP$spec) else spec
  )
  if(is.null(init)){
   parameters <- list(
     accuracy_coef   = rep(0,ngroup),
     threshold_coef  = rep(0,ngroup),
     shape_coef      = rep(0,ngroup),
     log_sigma_alpha = 0,
     log_sigma_theta = 0,
     alpha = rep(0, n_study),
     theta = rep(0, n_study)
    )
  } else {
    parameters <- init
  }
  
  
  dat2$model <- "RutterGatsonisReg"

  
  obj <- TMB::MakeADFun(dat2,
                        parameters,
                        map=if(is.null(map)) NULL else map,
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
  rep2 <- summary(rep,select="report")#p.value=TRUE)
  
  lspec <- length(dat2$spec)
  ## Get sensitivities and specificities
  qq   <- stats::qnorm(1-(1-conflevel)/2)
  rlse <- which(rownames(rep2)=="logitsens")
  sesp <- data.frame(spec=rep(dat2$spec,nrow(Z_pred)),
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

  res <- list(
    data         = XP,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    sensspec     = sesp  
  )
  class(res) <- c("RutterGatsonisReg")
  return(res)
}

