#' Fit the Nyaga Arm-Based DTA Network Meta-Analysis Model
#'
#' Fits the reduced arm-based ANOVA model for network meta-analysis of
#' diagnostic test accuracy proposed by Nyaga et al. The model uses a
#' binomial-normal likelihood and is fitted by maximum likelihood using a
#' custom \code{TMB} template and the Laplace approximation.
#'
#' The model permits multiple diagnostic tests to be evaluated within the
#' same study. Test arms from the same study share correlated study-level
#' random effects for logit-sensitivity and logit-specificity. Additional
#' independent arm-level random effects account for residual heterogeneity
#' among study-test combinations.
#'
#' For study \eqn{i} and test \eqn{k}, the model is
#'
#' \deqn{
#' \mathrm{logit}(Se_{ik})
#' =
#' \mu_{A,k}+u_{A,i}+\delta_{A,ik}
#' }
#'
#' and
#'
#' \deqn{
#' \mathrm{logit}(Sp_{ik})
#' =
#' \mu_{B,k}+u_{B,i}+\delta_{B,ik}.
#' }
#'
#' The study-level random effects follow
#'
#' \deqn{
#' (u_{A,i},u_{B,i})^\mathsf{T}
#' \sim
#' N(0,\Sigma),
#' }
#'
#' where \eqn{\Sigma} is a common unstructured variance-covariance matrix.
#' The arm-level effects are assumed independent, with
#'
#' \deqn{
#' \delta_{A,ik}\sim N(0,\tau_A^2)
#' }
#'
#' and
#'
#' \deqn{
#' \delta_{B,ik}\sim N(0,\tau_B^2).
#' }
#'
#' The summary sensitivity and specificity for test \eqn{k} are defined as
#' \eqn{\mathrm{logit}^{-1}(\mu_{A,k})} and
#' \eqn{\mathrm{logit}^{-1}(\mu_{B,k})}, respectively.
#'
#' @param data A \code{data.frame} containing one row for each observed
#'   study-test arm. Each study-test combination must occur only once.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name). Multiple rows may have the
#'   same study identifier when a study evaluates more than one diagnostic
#'   test.
#' @param test Diagnostic test identifier (column name). At least two
#'   diagnostic tests are required.
#' @param constrain Optional character string specifying a simplified
#'   covariance structure for the shared study-level random effects.
#'
#'   This can be useful for sparse networks, models with near-boundary
#'   variance estimates, or sensitivity analyses using reduced
#'   study-level heterogeneity structures.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       The standard unstructured study-level variance-covariance matrix
#'       is estimated. This is the default.
#'     }
#'     \item{\code{"sigma_AB"}}{
#'       The study-level covariance between logit-sensitivity and
#'       logit-specificity random effects is fixed at zero. The study-level
#'       sensitivity and specificity random effects remain independent.
#'     }
#'     \item{\code{"sigma2_A"}}{
#'       The study-level variance of logit-sensitivity is fixed at zero.
#'       This also implies a zero study-level covariance.
#'     }
#'     \item{\code{"sigma2_B"}}{
#'       The study-level variance of logit-specificity is fixed at zero.
#'       This also implies a zero study-level covariance.
#'     }
#'     \item{\code{"all"}}{
#'       Both study-level variances and their covariance are fixed at zero.
#'       Arm-level heterogeneity remains in the model unless it is
#'       separately removed using \code{tau_constrain}.
#'     }
#'   }
#' @param tau_constrain Optional character vector specifying arm-level
#'   variance components to be fixed at zero.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       Both arm-level variance components, \eqn{\tau_A^2} and
#'       \eqn{\tau_B^2}, are estimated. This is the default.
#'     }
#'     \item{\code{"tau2_A"}}{
#'       The arm-level variance of logit-sensitivity,
#'       \eqn{\tau_A^2}, is fixed at zero.
#'     }
#'     \item{\code{"tau2_B"}}{
#'       The arm-level variance of logit-specificity,
#'       \eqn{\tau_B^2}, is fixed at zero.
#'     }
#'   }
#'
#'   Both arm-level variance components may be fixed at zero by specifying
#'   \code{tau_constrain = c("tau2_A", "tau2_B")}. In that case, all
#'   residual heterogeneity among study-test arms is removed, while the
#'   shared study-level variance-covariance structure remains unless it is
#'   also restricted using \code{constrain}.
#' @param sensspec_constrain Optional character vector specifying
#'   restrictions on test-specific logit-sensitivity and/or
#'   logit-specificity parameters.
#'
#'   This can be used to fit reduced models in which summary sensitivity,
#'   summary specificity, or both are assumed equal across diagnostic tests.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       Separate logit-sensitivity and logit-specificity parameters are
#'       estimated for every diagnostic test. This is the default.
#'     }
#'     \item{\code{"sens"}}{
#'       All test-specific logit-sensitivity parameters are constrained to
#'       be equal. Test differences are therefore estimated only for
#'       logit-specificity.
#'     }
#'     \item{\code{"spec"}}{
#'       All test-specific logit-specificity parameters are constrained to
#'       be equal. Test differences are therefore estimated only for
#'       logit-sensitivity.
#'     }
#'   }
#'
#'   Both constraints may be specified simultaneously, for example
#'   \code{sensspec_constrain = c("sens", "spec")}. This constrains both
#'   summary sensitivity and summary specificity to be equal across tests.
#' @param conflevel Confidence level for confidence intervals. The default
#'   is \code{0.95}.
#' @param verbose Whether output from the inner
#'   \code{TMB} random-effects optimization should be printed. The default
#'   is \code{FALSE}.
#'
#' @details
#' The input data must be in arm-based long format, with one row for each
#' observed combination of study and diagnostic test. Studies evaluating
#' multiple tests should therefore contribute multiple rows with the same
#' study identifier.
#'
#' The model estimates one common study-level variance-covariance matrix
#' \eqn{\Sigma} across all diagnostic tests. Separate arm-level variance
#' components \eqn{\tau_A^2} and \eqn{\tau_B^2} describe residual
#' heterogeneity in logit-sensitivity and logit-specificity among
#' study-test arms.
#'
#' Estimates of a study-level or arm-level variance component may occur at
#' or near zero. In that situation, Wald standard errors for that variance
#' component and related quantities may be unreliable, although inference
#' for test-specific fixed effects may remain usable. The constraint
#' arguments allow reduced models to be fitted explicitly when supported by
#' the data or required for sensitivity analyses.
#'
#' @return A list of class \code{"Nyaga"} and \code{"DTAmodel"} with
#' components:
#' \itemize{
#'   \item \code{data}: The analyzed data with derived quantities.
#'   \item \code{fit}: Optimization result returned by
#'   \code{\link[stats]{nlminb}}.
#'   \item \code{sdreport}: The \code{TMB} standard report.
#'   \item \code{sdreport2}: Summary of reported and derived model
#'   parameters with test-specific row names.
#'   \item \code{vcov}: Variance-covariance matrix of reported parameters.
#'   \item \code{sensspec}: Test-specific conditional sensitivity and
#'   specificity estimates with confidence intervals.
#'   \item \code{LRDOR}: Test-specific positive and negative likelihood
#'   ratios and diagnostic odds ratios with confidence intervals.
#'   \item \code{tests}: Diagnostic test levels used in the fitted model,
#'   in model-matrix order. The first level is the reference test.
#'   \item \code{constrain}: Study-level random-effects parameters fixed
#'   at zero.
#'   \item \code{tau_constrain}: Arm-level variance components fixed at
#'   zero.
#'   \item \code{sensspec_constrain}: Whether sensitivity and/or
#'   specificity were constrained to be common across tests.
#' }
#'
#' @examples
#' \dontrun{
#' data("HPV")
#'
#' fit <- fitNyaga(
#'   data = HPV,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = StudyID,
#'   test = Test
#' )
#'
#' fit
#' summary(fit)
#' }
#'
#' @references
#' Nyaga, V. N., Aerts, M., & Arbyn, M. (2018).
#' ANOVA model for network meta-analysis of diagnostic test accuracy data.
#' \emph{Statistical Methods in Medical Research}, 27(6), 1766--1784.
#' \doi{10.1177/0962280216669182}
#'
#' Nyaga, V. N., & Arbyn, M. (2023).
#' Comparison and validation of \emph{metadta} for meta-analysis of
#' diagnostic test accuracy studies.
#' \emph{Research Synthesis Methods}, 14(3), 544--562.
#' \doi{10.1002/jrsm.1634}
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb qnorm plogis
#' @export

fitNyaga <- function(data,
                     TP, 
                     FP, 
                     FN, 
                     TN,
                     study,
                     test,
                     constrain=NULL,
                     tau_constrain=NULL,
                     sensspec_constrain=NULL,
                     conflevel=0.95,
                     verbose=FALSE) {
  
  allowed_sensspec_constraints <- c(
    "sens",
    "spec"
  )
  if (length(sensspec_constrain) == 0L) {
    sensspec_constrain <- NULL
  }
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
  
  ## Constraints ###
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
  
  allowed_tau_constraints <- c(
    "tau2_A",
    "tau2_B"
  )
  
  if (length(tau_constrain) == 0L) {
    tau_constrain <- NULL
  }
  
  if (!is.null(tau_constrain)) {
    if (!is.character(tau_constrain)) {
      stop("'tau_constrain' must be a character vector or NULL."
      )
    }
    
    tau_constrain <- unique(tau_constrain)
    
    invalid_tau_constraints <- setdiff(
      tau_constrain,
      allowed_tau_constraints
    )
    
    if (length(invalid_tau_constraints) > 0L) {
      stop(
        "Unknown arm-level variance constraint(s): ",
        paste(invalid_tau_constraints, collapse = ", "),
        ". Allowed values are 'tau2_A' and 'tau2_B'."
      )
    }
  }
  
  TP_col       <- deparse(substitute(TP))
  FP_col       <- deparse(substitute(FP))
  FN_col       <- deparse(substitute(FN))
  TN_col       <- deparse(substitute(TN))
  study_col    <- deparse(substitute(study))
  test_col     <- deparse(substitute(test))
  
  X <- data.frame(study = data[[study_col]],
                  TP = data[[TP_col]],
                  FP = data[[FP_col]],
                  FN = data[[FN_col]],
                  TN = data[[TN_col]],
                  test = data[[test_col]])
  
  X <- check_data(dat=X,conflevel=conflevel)
  
  X$test  <- droplevels(factor(X$test))
  X$study <- droplevels(factor(X$study))
  ltest   <- levels(X$test)
  n_arm   <- nrow(X)
  n_study <- nlevels(X$study)
  n_test  <- length(ltest)
  if (n_test < 2L) { stop( "'fitNyaga()' requires at least two diagnostic tests." ) }
  if (anyDuplicated(X[c("study", "test")])) {
    stop(
      "Each study-test combination must occur only once. ",
      "Aggregate duplicate rows to one 2x2 table per study-test arm."
    )
  }
  
  XP <- getXP(X=X)
  
  ### Get initial values
  init <- tryCatch(suppressWarnings(
                   fitReitsmaSubgroupTMB(data=X,
                                         TP=TP,
                                         FP=FP,
                                         FN=FN,
                                         TN=TN,
                                         study=study,
                                         subgroup=test,
                                         conflevel=conflevel,
                                         variances="common",
                                         constrain=NULL)$sdreport),
                   error = function(e) NULL)
  
  muA_init <- rep(NA_real_,n_test)
  muB_init <- rep(NA_real_,n_test)
  lsA_init <- NA_real_
  lsB_init <- NA_real_
  thetaAB_init <- NA_real_
  ltA_init <- NA_real_
  ltB_init <- NA_real_
  
  sigma_fraction <- 0.7
  
  if(!is.null(init)){
    muA_init      <- init$par.fixed[names(init$par.fixed) == "mu_A"]
    muB_init      <- init$par.fixed[names(init$par.fixed) == "mu_B"]
    lsA_init      <- init$par.fixed["log_sigma_A"] + 0.5*log(sigma_fraction)
    lsB_init      <- init$par.fixed["log_sigma_B"] + 0.5*log(sigma_fraction)
    thetaAB_init  <- init$par.fixed["theta_AB"]
    ltA_init      <- init$par.fixed["log_sigma_A"] + 0.5*log(1-sigma_fraction)
    ltB_init      <- init$par.fixed["log_sigma_B"] + 0.5*log(1-sigma_fraction)
  }
  
  muA_init[!is.finite(muA_init)] <- 0
  muB_init[!is.finite(muB_init)] <- 0
  if (!is.finite(lsA_init)) {lsA_init <- log(0.5)}
  if (!is.finite(lsB_init)) {lsB_init <- log(0.5)}
  if (!is.finite(thetaAB_init)) {thetaAB_init <- 0}
  if (!is.finite(ltA_init)) {ltA_init <- log(0.25)}
  if (!is.finite(ltB_init)) {ltB_init <- log(0.25)}
  
  parameters <- list(
    mu_A = muA_init,
    mu_B = muB_init,
    
    log_sigma_A = lsA_init,
    log_sigma_B = lsB_init,
    theta_AB = thetaAB_init,
    
    log_tau_A = ltA_init,
    log_tau_B = ltB_init,
    
    sensu = rep(0, n_study),
    specu = rep(0, n_study),
    
    sens_delta = rep(0, n_arm),
    spec_delta = rep(0, n_arm)
  )
  
  dat2 <- list(
    TP    = X$TP,
    FP    = X$FP,
    FN    = X$FN,
    TN    = X$TN,
    study = as.numeric(X$study) - 1,
    test  = as.numeric(X$test) - 1
  )
  
  ## Constraints ###
  map <- list()
  
  if (!is.null(sensspec_constrain)) {
    if ("sens" %in% sensspec_constrain) {
      map$mu_A <- factor(rep(1,n_test))
    }
    if ("spec" %in% sensspec_constrain) {
      map$mu_B <- factor(rep(1,n_test))
    }
  }
  
  if(!is.null(constrain)){
    if (constrain == "sigma_AB") {
      parameters$theta_AB <- 0
      map$theta_AB        <- factor(NA)
    }
    if (constrain == "sigma2_A") {
      parameters$log_sigma_A <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_A        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
    if (constrain == "sigma2_B") {
      parameters$log_sigma_B <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_B        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
    if (constrain == "all") {
      parameters$log_sigma_A <- log(.Machine$double.eps)
      parameters$log_sigma_B <- log(.Machine$double.eps)
      parameters$theta_AB    <- 0
      map$log_sigma_A        <- factor(NA)
      map$log_sigma_B        <- factor(NA)
      map$theta_AB           <- factor(NA)
    }
  }
  
  if (!is.null(tau_constrain)) {
    if ("tau2_A" %in% tau_constrain) {
      parameters$log_tau_A <- log(.Machine$double.eps)
      map$log_tau_A        <- factor(NA)
    }
    if ("tau2_B" %in% tau_constrain) {
      parameters$log_tau_B <- log(.Machine$double.eps)
      map$log_tau_B        <- factor(NA)
    }
  }
  
  
  ###
  
  dat2$model = "Nyaga"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters=parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("sensu","specu","sens_delta","spec_delta"),
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
  if(!rep$pdHess){
    warning(
      "The Hessian matrix is not positive definite. ",
      "Parameter estimates and standard errors may be unreliable. ",
      "Consider checking model specification, data quality, or potential identifiability issues."
    )
  }
  rep1 <- summary(rep, select = "report")
  rep2 <- rename_nyaga_rows(rep1,ltest)
  
  
  # Variance covariance matrix of fixed effects
  vcov <- rep$cov
  colnames(vcov) <- rownames(vcov) <- rownames(rep2)
  
  ### Sensitivity and Specificity
  qq             <- stats::qnorm(1-(1-conflevel)/2)
  sesp           <- as.data.frame(rep2[grepl("^mu_",rownames(rep2)),,drop=FALSE])
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,c("Orig","conflevel","CI_Lower","CI_Upper")]
  colnames(sesp) <- c("Estimate","conflevel","CI_Lower","CI_Upper")
  sesp$type      <- sub("^mu_([AB])\\..*$", "\\1", rownames(sesp))
  sesp$type      <- c(A = "sens", B = "spec")[sesp$type]
  sesp           <- sesp[,c("type","Estimate","conflevel","CI_Lower","CI_Upper")]
  ###
  lrdor2 <- data.frame()
  for(i in seq_along(ltest)){
    sg      <- ltest[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    lsens   <- rep2[mu_A.sg,"Estimate"]
    lspec   <- rep2[mu_B.sg,"Estimate"]
    S       <- vcov[c(mu_A.sg,mu_B.sg),c(mu_A.sg,mu_B.sg)]
    lrdor   <- getLRDOR(lsens=lsens, lspec=lspec, S=S, conflevel=conflevel)
    rownames(lrdor) <- paste0(ltest[i],": ",rownames(lrdor))
    lrdor2  <- rbind(lrdor2,lrdor)
  }
  
  ##
  ret <- list(data         = XP,
              fit          = fit,
              sdreport     = rep,
              sdreport2    = rep2,
              vcov         = vcov,
              sensspec     = sesp,
              LRDOR        = lrdor2,
              tests        = ltest,
              constrain    = constrain,
              tau_constrain=tau_constrain,
              sensspec_constrain=sensspec_constrain)
  class(ret) <- c("Nyaga","DTAmodel")
  return(ret)
}


#' @keywords internal
#' @noRd
rename_nyaga_rows <- function(ss, ltest) {
  rn <- rownames(ss)
  if (is.null(rn)) rn <- rep("", nrow(ss))
  
  i_nuA   <- which(rn == "nu_A")
  i_nuB   <- which(rn == "nu_B")
  
  nu_A_names  <- paste0("nu_A.", ltest[-1])
  nu_B_names  <- paste0("nu_B.", ltest[-1])
  
  i_muA   <- which(rn == "mu_A")
  i_muB   <- which(rn == "mu_B")
  
  mu_A_names  <- paste0("mu_A.", ltest)
  mu_B_names  <- paste0("mu_B.", ltest)
  
  rn[i_nuA] <- nu_A_names
  rn[i_nuB] <- nu_B_names
  
  rn[i_muA] <- mu_A_names
  rn[i_muB] <- mu_B_names
  
  rownames(ss) <- rn
  ss
}


#' @rdname logLik.dtametaTMB
#' @export
logLik.Nyaga <- logLik.DTAmodel

#' @rdname coef.dtametaTMB
#' @export
coef.Nyaga <- coef.DTAmodel

#' @rdname vcov.dtametaTMB
#' @export
vcov.Nyaga <- vcov.DTAmodel


#' Print Method for Nyaga Objects
#'
#' Displays a concise summary of a fitted Nyaga arm-based diagnostic test
#' accuracy network meta-analysis model, including number of studies, convergence status
#' and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"Nyaga"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.Nyaga}}
#' @return
#' Invisibly returns the input object.
#' @method print Nyaga
#' @export
print.Nyaga <- function(x, ...) {
  
  cat("\n", "Nyaga NMA Model", "\n", sep = "")
  cat(strrep("-", nchar("Nyaga NMA Model")), "\n\n", sep = "")
  
  n_study   <- nlevels(x$data$study)
  n_arm     <- nrow(x$data)
  n_test    <- length(x$tests)
  converged <- x$fit$convergence == 0
  pdHess    <- x$sdreport$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies   :", n_study, "\n")
  cat("Number of arms      :", n_arm, "\n")
  cat("Number of tests     :", n_test, "\n")
  cat("Optimizer           :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian             :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|          :", max(abs(x$sdreport$gradient.fixed)), "\n")
  cat("-2 log likelihood   :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC                 :", round(AIC(x), 3), "\n")
  cat("BIC                 :", round(BIC(x), 3), "\n")
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}


#' Summary Method for Nyaga Objects
#'
#' Extracts key results from an object of class \code{"Nyaga"}, including
#' parameter estimates, sensitivity and specificity summaries.
#'
#' @param object An object of class \code{"Nyaga"} as returned by
#'   \code{\link{fitNyaga}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class \code{"summary.Nyaga"} with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{LRDOR}: Estimated diagnostic odds ratio and likelihood ratios with confidence intervals.
#'   \item \code{tests}: Test names.
#' }
#'
#' @seealso \code{\link{fitNyaga}}
#' @method summary Nyaga
#' @export
summary.Nyaga <- function(object, ...) {
  ret <- list(estimates = object$sdreport2,
              sensspec = object$sensspec,
              LRDOR = object$LRDOR,
              tests = object$tests)
  class(ret) <- c("summary.Nyaga","summary.DTAmodel")
  return(ret)
}




#' Plot Results from a Nyaga Model
#'
#' Produces a summary ROC plot for objects of class \code{"Nyaga"}
#' obtained from \code{\link{fitNyaga}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. 
#'
#' @param x An object of class \code{"Nyaga"}, as returned by
#'   \code{\link{fitNyaga}}.
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the tests. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally. Relative adjustments can be 
#' specified directly, for example, \code{scale = 0.02 * 0.3} produces marker 
#' dimensions equal to 30 percent of the default.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"equal"}{All studies shown with equal size. Default.}
#'    \item{"sampsize"}{
#'     Horizontal and vertical symbol dimensions are proportional to the
#'     relative numbers of non-diseased and diseased participants,
#'     respectively.}
#'    \item{"se"}{
#'       Horizontal and vertical symbol dimensions are proportional to
#'       normalized approximate inverse standard errors of observed logit
#'       specificity and logit sensitivity, respectively.}
#'    \item{"sampsize_revman"}{
#'       RevMan-style sample-size scaling. Horizontal and vertical symbol
#'       dimensions are based on the square roots of the numbers of
#'       non-diseased and diseased participants, respectively. A common
#'       analysis-specific scaling constant is chosen so that the average
#'       larger dimension measure across studies equals nine before adding
#'       the fixed baseline dimension.}
#'    \item{"se_revman"}{
#'       RevMan-style inverse-standard-error scaling. Horizontal and
#'       vertical symbol dimensions are based on approximate inverse standard
#'       errors of observed logit specificity and logit sensitivity,
#'       respectively. A common analysis-specific scaling constant is chosen
#'       so that the average larger dimension measure across studies equals
#'       nine before adding the fixed baseline dimension.}
#'  } 
#' @param col Vector of colours used for test-specific
#'   study-level symbols, summary points, confidence and prediction region. 
#'   If \code{NULL}, colours are generated automatically.
#' @param nudge_legend Numeric horizontal offset for the test legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Network Meta-Analysis"}.
#' 
#' @param conflevel Confidence level for the confidence region. Default is \code{0.95}.
#' @param predlevel Confidence level for the prediction region. Default is \code{0.95}.
#' @param connectstudies Whether the point estimates of two diagnostic
#'   tests evaluated within the same study should be connected.
#'   
#' @param studyCI Whether study-level confidence
#'   intervals should be displayed for the observed sensitivity and
#'   specificity estimates. Defaults to \code{FALSE}.
#'
#' @param studylabels Whether study identifiers should
#'   be displayed next to the study-level estimates. Defaults to
#'   \code{FALSE}.
#'   
#' @param legend Whether a legend should be displayed. Default is \code{TRUE}.
#'
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The plot is constructed on the ROC scale with sensitivity on the y-axis
#' and specificity on the x-axis (displayed as 1 - false positive rate on a reversed axis).
#'
#' Study-specific estimates are shown with study-specific symbols and colours.
#'
#' The following elements are displayed:
#' \itemize{
#'   \item Study-level sensitivity and specificity estimates
#'   \item Summary (pooled) estimate
#'   \item confidence region around the summary point
#'   \item prediction region reflecting between-study variability
#' }
#' 
#' Confidence and prediction contours use the finite-sample F-based convention 
#' employed by other plotting methods in dtametaTMB. Prediction regions incorporate 
#' both shared study-level and arm-level heterogeneity.
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitNyaga}}
#' @method plot Nyaga
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.Nyaga <- function(x,
                       symbol=NULL,
                       scale=0.02, 
                       size=c("equal","sampsize","se","sampsize_revman","se_revman"), 
                       main="Diagnostic Test Accuracy Network Meta-Analysis",
                       col=NULL,
                       nudge_legend=-0.4,
                       conflevel=0.95,
                       predlevel=0.95,
                       connectstudies=FALSE,
                       studyCI=FALSE,
                       studylabels=FALSE,
                       legend=TRUE,
                       ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  if (!is.numeric(predlevel) || length(predlevel) != 1L ||
      predlevel <= 0 || predlevel >= 1) {
    stop("predlevel must be a single number in (0, 1).")
  }
  if(connectstudies) {
    if(length(unique(x$data$test)) != 2) {
      warning("'connectstudies=TRUE' is only recommended for two-test comparisons." )
    }
  }
  size   <- match.arg(size)
  test   <- levels(x$data$test)
  ntest  <- length(test)
  warn_unestimable_sroc_points(x$data)
  aes  <- validateSubgroupAesthetics(col=col,symbol=symbol,nsub=ntest)
  col  <- aes$col
  col2 <- aes$col2
  symbols2 <- aes$symbols2
  symb <- vector(mode="integer",length=ntest)
  # Calculations for percentage weights
  pct <- getWEIGHTS(xdata=x$data,size=size)
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  for (i in seq_along(test)){
    sg    <- x$data$test==test[i]
    xspec <- x$data$spec[sg]
    xsens <- x$data$sens[sg]
    pctsp <- pct$sp[sg]
    pctse <- pct$se[sg]
    xstud <- x$data$study[sg]
    xdata <- x$data[sg,]
    ####
    symb[i] <- plotStudyLevelEstimates(xsens=xsens,xspec=xspec,xstud=xstud,xdata=xdata,
                                       symbol=symbols2[i],pctse=pctse,pctsp=pctsp,
                                       scale=scale,col=col2[i],
                                       conflevel=conflevel,studyCI=studyCI,
                                       studylabels=studylabels,LCA=FALSE)
  }
  # Add summary point
  for (i in seq_along(test)){
    sg      <- test[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    mean_point <- data.frame(1-x$sensspec[mu_B.sg,"Estimate"],
                             x$sensspec[mu_A.sg,"Estimate"])
    points(mean_point, col=col[i], cex=1.5, pch=15)
  }
  # Add confidence and prediction region
  for(i in seq_along(test)){
    sg      <- test[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg) 
    s2_A.sg <- "total_var_A"
    s2_B.sg <- "total_var_B"
    s_AB.sg <- "sigma_AB"
    muA     <- x$sdreport2[mu_A.sg,"Estimate"]
    muB     <- x$sdreport2[mu_B.sg,"Estimate"]
    seA     <- x$sdreport2[mu_A.sg,"Std. Error"]
    seB     <- x$sdreport2[mu_B.sg,"Std. Error"]
    covAB   <- x$vcov[mu_A.sg,mu_B.sg]
    varA    <- x$sdreport2[s2_A.sg,"Estimate"]
    varB    <- x$sdreport2[s2_B.sg,"Estimate"]
    sAB     <- x$sdreport2[s_AB.sg,"Estimate"]
    region  <- getConfPredRegion(muA=muA,muB=muB,
                                 seA=seA,seB=seB,covAB=covAB, # conf
                                 varA=varA,varB=varB,sAB=sAB, # pred
                                 nstudy=sum(x$data$test == sg, na.rm = TRUE),
                                 conflevel=conflevel,
                                 predlevel=predlevel)
    lines(region$conf, lty=2, lwd=2, col=col2[i])
    lines(region$pred, lty=3, lwd=2, col=col2[i])
  }
  ## Connect studies
  if(isTRUE(connectstudies)){
    for(st in unique(x$data$study)){
      tmp <- x$data[x$data$study == st, ]
      tmp <- tmp[order(tmp$test), ]
      if(nrow(tmp) == 2) {
        lines(x = c(1-tmp$spec[1],1-tmp$spec[2]),
              y = c(tmp$sens[1],tmp$sens[2]),
              type="l",col = "grey70",lwd =1)
      }
    }
  }
  # Add the legend 
  conf_lab <- paste0(round(100 * conflevel), "% Confidence region")
  pred_lab <- paste0(round(100 * predlevel), "% Prediction region")
  if(isTRUE(legend)){
      graphics::legend("bottomright", 
                       bty ="n",
                       legend = c(NA,
                                  "Summary estimate",
                                  conf_lab,
                                  pred_lab,
                                  "Study estimates"), 
                       pch = c(NA,15,NA,NA,symb[1]), 
                       lty = c(NA,NA,2,3,NA), 
                       lwd = c(NA,NA,2,2,NA), 
                       col = c(NA,"black","black","black","darkgray"))
  }
  graphics::legend("right",
                   inset = c(nudge_legend, 0),
                   legend = test,
                   col    = col,
                   pch = symb,
                   xpd = TRUE,
                   cex = 1.2,
                   bty = "n")
  invisible(NULL)
}



#' Coupled Forest Plot for Diagnostic Test Accuracy Meta-Analysis
#' 
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#'
#' @param x Object of class \code{"Nyaga"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param test_label Column name for the test. Defaults to \code{"Test"}.
#' @param order Specifies the ordering of studies
#'   in the forest plot. Can be `"study"` (default), which orders
#'   observations by study identifier and then subgroup, or
#'   `"test"`, which orders observations by test and then
#'   study identifier.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest Nyaga
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.Nyaga <- function(x, conflevel=0.95, test_label="Test", order=c("study","test"), ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  orderby <- match.arg(order)
  ss <- getForestSensSpec(xdata=x$data,conflevel=conflevel)
  
  XP <- ss$XP
  if(orderby=="study"){XP <- XP[base::order(XP$study,XP$test), ]}
  if(orderby=="test"){XP <- XP[base::order(XP$test,XP$study), ]}
  dt <- XP[,c("study","test","TP","FP","FN","TN","senslabel","speclabel")]
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  cc <- colnames(dt) 
  colnames(dt) <- c("Study",test_label,cc[3:6],
                    ss$senslab,ss$speclab," ",ss$senslab," ",ss$speclab)
  
  getForestPlot(dt=dt,XP=XP,r=1,s=3)
  
}

