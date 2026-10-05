#' Fit the Reitsma Subgroup TMB Model
#'
#' Fits the Reitsma bivariate random-effects model with a single categorical covariate
#' for diagnostic test accuracy (DTA) meta-analysis using a binomial-normal 
#' likelihood via a custom \code{TMB} template.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param subgroup A single categorical study-level subgroup variable (column name).
#' @param constrain Optional character string specifying a simplified
#'   covariance structure for the Reitsma model.
#'
#'   This can be useful for sparse data, small meta-analyses, or for
#'   reproducing simplified bivariate models described in the Cochrane
#'   Handbook for Diagnostic Test Accuracy Reviews.
#'
#'   Allowed values are:
#'   \describe{
#'     \item{\code{NULL}}{
#'       The standard unconstrained Reitsma model is fitted.
#'     }
#'     \item{\code{"sigma_AB"}}{
#'       The covariance between logit-sensitivity and
#'       logit-specificity random effects is fixed at zero.
#'       Random effects remain independent.
#'     }
#'     \item{\code{"sigma2_A"}}{
#'       The between-study variance of logit-sensitivity is fixed at zero.
#'       This also implies a zero covariance.
#'     }
#'     \item{\code{"sigma2_B"}}{
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
#'   logit-specificity parameters.
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
#'   matrix should be assumed to be \code{"common"} (default) or \code{"unequal"} 
#'   across subgroups. If \code{"common"}, a single between-study variance-covariance 
#'   matrix is estimated and shared across all subgroups. If \code{"unequal"}, 
#'   subgroup-specific variance-covariance matrices are estimated.
#'   
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
#'
#' @return A list of class \code{"ReitsmaSubgroupTMB"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities
#'   \item \code{fit}: Optimization result from \code{nlminb}.
#'   \item \code{sdreport}: TMB standard report.
#'   \item \code{sdreport2}: Summary of reported parameters.
#'   \item \code{vcov}: variance-covariance matrix of reported parameters.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{LRDOR}: Diagnostic odds ratios and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{subgroups}: The subgroup levels used in the model fit.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#'   \item \code{sensspec_constrain}: Whether sensitivities/specificities are assumed common across subgroups.
#'   \item \code{variances}: Variance structure used in the fitted model.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsmaSubgroupTMB(
#'   data = anticcp,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study,
#'   subgroup = generation
#' )
#' fit
#' summary(fit)
#' 
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
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats binomial complete.cases qnorm plogis qlogis vcov cor sd reshape
#' @export
fitReitsmaSubgroupTMB <- function(data,
                                  TP, FP, FN, TN,
                                  study,
                                  subgroup,
                                  constrain=NULL,
                                  sensspec_constrain=NULL,
                                  variances=c("common","unequal"),
                                  conflevel=0.95,
                                  verbose=FALSE) {
  
  variances  <- match.arg(variances)
  
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
  
  X <- check_data(dat=X,
                  conflevel=conflevel)
  
  X$subgroup <- droplevels(factor(X$subgroup))
  lsub       <- levels(X$subgroup)
  n_study    <- nrow(X)
  G          <- length(lsub)
  
  XP <- getXP(X=X)
  
  ### Get initial values
  init <- fitReitsmaTMB(data=X,
                        TP=TP,FP=FP,FN=FN,TN=TN,study=study,
                        constrain=NULL,
                        conflevel=conflevel)$sdreport2
  muA_init     <- init["mu_A.sens","Estimate"]
  muB_init     <- init["mu_B.spec","Estimate"]
  sA_init      <- sqrt(init["sigma2_A.sens","Estimate"])
  sA_init      <- max(sA_init,1e-05)
  sB_init      <- sqrt(init["sigma2_B.spec","Estimate"])
  sB_init      <- max(sB_init,1e-05)
  sAB_init     <- init["sigma_AB","Estimate"]
  rAB_init     <- sAB_init/(sA_init*sB_init)
  rAB_init     <- max(min(rAB_init,0.99),-0.99)
  if(is.na(rAB_init)) rAB_init <- 0

  ###
  parameters <- list(
    mu_A = rep(muA_init,G),
    mu_B = rep(muB_init,G),
    
    log_sigma_A  = rep(log(sA_init),G),
    log_sigma_B  = rep(log(sB_init),G),
    theta_AB     = rep(atanh(rAB_init),G),
    
    sensu = rep(0,n_study),
    specu = rep(0,n_study)
  )
  
  dat2 <- list(
    TP = X$TP,
    FP = X$FP,
    FN = X$FN,
    TN = X$TN,
    group = as.numeric(X$subgroup)-1
  )
  
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
  
  map <- list()
  
  if(variances=="common"){
    map$log_sigma_A <- factor(rep(1,G))
    map$log_sigma_B <- factor(rep(1,G))
    map$theta_AB    <- factor(rep(1,G))
  }
  
  if(variances=="common"){
    if(!is.null(constrain)){
      if (constrain == "sigma_AB") {
        parameters$theta_AB    <- rep(0,G)
        map$theta_AB           <- rep(factor(NA),G)
      }
      if (constrain == "sigma2_A") {
        parameters$log_sigma_A <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB    <- rep(0,G)
        map$log_sigma_A        <- rep(factor(NA),G)
        map$theta_AB           <- rep(factor(NA),G)
      }
      if (constrain == "sigma2_B") {
        parameters$log_sigma_B <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB    <- rep(0,G)
        map$log_sigma_B        <- rep(factor(NA),G)
        map$theta_AB           <- rep(factor(NA),G)
      }
      if (constrain == "all") {
        parameters$log_sigma_A <- rep(log(.Machine$double.eps),G)
        parameters$log_sigma_B <- rep(log(.Machine$double.eps),G)
        parameters$theta_AB    <- rep(0,G)
        map$log_sigma_A        <- rep(factor(NA),G)
        map$log_sigma_B        <- rep(factor(NA),G)
        map$theta_AB           <- rep(factor(NA),G)
      }
    }
  }
  
  if (!is.null(sensspec_constrain)) {
    if ("sens" %in% sensspec_constrain) {
      map$mu_A <- factor(rep(1,G))
    }
    if ("spec" %in% sensspec_constrain) {
      map$mu_B <- factor(rep(1,G))
    }
  }
  
  dat2$model = "ReitsmaSubgroup"
  
  # TMB Objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        map = if(length(map) == 0) NULL else map,
                        random = c("sensu", "specu"),
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
  rep2 <- rename_reitsub_rows(rep1,lsub)

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
  ### Diagnostic odds ratios and Likelihood ratios
  lrdor2 <- data.frame()
  for(i in seq_along(lsub)){
    sg      <- lsub[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
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
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    s2_A.sg <- paste0("sigma2_A.",sg)
    s2_B.sg <- paste0("sigma2_B.",sg)
    s_AB.sg <- paste0("sigma_AB.",sg)
    ruga <- getRUGA(lsens=rep2[mu_A.sg,"Estimate"],
                    lspec=rep2[mu_B.sg,"Estimate"],
                    sigma_a=sqrt(rep2[s2_A.sg,"Estimate"]),
                    sigma_b=sqrt(rep2[s2_B.sg,"Estimate"]),
                    sigma_ab=rep2[s_AB.sg,"Estimate"])
    ruga2 <- rbind(ruga2,ruga)
  }
  rownames(ruga2) <- lsub
  
  ##
  ret <- list(data         = XP,
              fit          = fit,
              sdreport     = rep,
              sdreport2    = rep2,
              vcov         = vcov,
              sensspec     = sesp,
              LRDOR        = lrdor2,
              RutterGatsonis_recovered = ruga2,
              subgroups    = lsub,
              constrain    = constrain,
              sensspec_constrain = sensspec_constrain,
              variances    = variances)
  class(ret) <- c("ReitsmaSubgroupTMB","Reitsma","DTAmodel")
  return(ret)
}

#' @keywords internal
#' @noRd
rename_reitsub_rows <- function(ss, lsub) {
  
  rn <- rownames(ss)
  if (is.null(rn)) rn <- rep("", nrow(ss))
  
  i_nuA   <- which(rn == "nu_A")
  i_nuB   <- which(rn == "nu_B")
  
  nu_A_names     <- paste0("nu_A.", lsub[-1])
  nu_B_names     <- paste0("nu_B.", lsub[-1])
  
  i_muA   <- which(rn == "mu_A")
  i_muB   <- which(rn == "mu_B")
  i_s2A   <- which(rn == "sigma2_A")
  i_s2B   <- which(rn == "sigma2_B")
  i_sAB   <- which(rn == "sigma_AB")
  i_rAB   <- which(rn == "rho_AB")
  
  mu_A_names        <- paste0("mu_A.", lsub)
  mu_B_names        <- paste0("mu_B.", lsub)
  sigma2_A_names    <- paste0("sigma2_A.", lsub)
  sigma2_B_names    <- paste0("sigma2_B.", lsub)
  sigma_AB_names    <- paste0("sigma_AB.", lsub)
  rho_AB_names      <- paste0("rho_AB.", lsub)
  
  rn[i_nuA] <- nu_A_names
  rn[i_nuB] <- nu_B_names
  
  rn[i_muA] <- mu_A_names
  rn[i_muB] <- mu_B_names
  rn[i_s2A] <- sigma2_A_names
  rn[i_s2B] <- sigma2_B_names
  rn[i_sAB] <- sigma_AB_names
  rn[i_rAB] <- rho_AB_names
  
  rownames(ss) <- rn
  ss
}



#' Print Method for ReitsmaSubgroupTMB Objects
#'
#' Displays a concise summary of a fitted ReitsmaSubgroupTMB diagnostic test
#' accuracy model, including number of studies, convergence status
#' and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"ReitsmaSubgroupTMB"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.ReitsmaSubgroupTMB}}
#' @return
#' Invisibly returns the input object.
#' @method print ReitsmaSubgroupTMB
#' @export
print.ReitsmaSubgroupTMB <- function(x, ...) {
  
  cat("\n", "Reitsma Subgroup TMB Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma Subgroup TMB Model")), "\n\n", sep = "")
  
  n_study   <- nrow(x$data)
  n_sub     <- length(x$subgroups)
  converged <- x$fit$convergence == 0
  pdHess    <- x$sdreport$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies   :", n_study, "\n")
  cat("Number of subgroups :", n_sub, "\n")
  cat("Optimizer           :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian             :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|          :", max(abs(x$sdreport$gradient.fixed)), "\n")
  cat("-2 log likelihood   :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC                 :", round(AIC(x), 3), "\n")
  cat("BIC                 :", round(BIC(x), 3), "\n")
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.ReitsmaSubgroupTMB <- logLik.DTAmodel

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaSubgroupTMB <- coef.DTAmodel


#' Coupled Forest Plot for Diagnostic Test Accuracy Meta-Analysis
#' 
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#'
#' @param x Object of class \code{"ReitsmaSubgroupTMB"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param subgroup_label Column name for the subgroup. Defaults to \code{"Subgroup"}.
#' @param order Specifies the ordering of studies
#'   in the forest plot. Can be `"study"` (default), which orders
#'   observations by study identifier and then subgroup, or
#'   `"subgroup"`, which orders observations by subgroup and then
#'   study identifier.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest ReitsmaSubgroupTMB
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qbeta
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.ReitsmaSubgroupTMB <- function(
    x,
    conflevel = 0.95,
    subgroup_label = "Subgroup",
    order = c("study", "subgroup"),
    ...) {
  
  forest.ReitsmaSubgroup(
    x = x,
    conflevel = conflevel,
    subgroup_label = subgroup_label,
    order = order,
    ...
  )
}



#' Summary Method for ReitsmaSubgroupTMB Objects
#'
#' Extracts key results from an object of class \code{"ReitsmaSubgroupTMB"}, including
#' parameter estimates, sensitivity and specificity summaries, and recovered
#' HSROC parameters from the Rutter–Gatsonis parameterization.
#'
#' @param object An object of class \code{"ReitsmaSubgroupTMB"} as returned by
#'   \code{\link{fitReitsmaSubgroupTMB}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class \code{"summary.ReitsmaSubgroupTMB"} with the following components:
#' \itemize{
#'   \item \code{estimates}: Parameter estimates with standard errors.
#'   \item \code{sensspec}: Estimated sensitivity and specificity with confidence intervals.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{subgroups} Subgroup names.
#' }
#'
#' @seealso \code{\link{fitReitsmaSubgroupTMB}}
#'
#'
#' @export
summary.ReitsmaSubgroupTMB <- function(object, ...) {
  ret <- list(estimates = object$sdreport2,
              sensspec = object$sensspec,
              RutterGatsonis_recovered = object$RutterGatsonis_recovered,
              subgroups = object$subgroups)
  class(ret) <- c("summary.ReitsmaSubgroupTMB","summary.Reitsma","summary.DTAmodel")
  return(ret)
}




#' Plot Results from a Reitsma Subgroup TMB Model
#'
#' Produces a summary ROC plot for objects of class \code{"ReitsmaSubgroupTMB"}
#' obtained from \code{\link{fitReitsmaSubgroupTMB}}. The plot displays study-level
#' estimates of sensitivity and specificity, the summary operating point,
#' and corresponding confidence and prediction regions. Optionally, the
#' HSROC (hierarchical summary ROC) curve can be overlaid.
#'
#' @param x An object of class \code{"ReitsmaSubgroupTMB"}, as returned by
#'   \code{\link{fitReitsmaSubgroupTMB}}.
#' @param symbol Character vector specifying the symbols used for study-level
#'   estimates of the subgroups. One can choose \code{"rectangle"}, \code{"ellipse"},
#'   \code{"diamond"}, \code{"triangle"}, \code{"cross"}, \code{"plus"},
#'   or \code{"star"}. If \code{NULL}, symbols are generated automatically.
#' @param scale Numeric scaling factor controlling study-symbol dimensions.
#' The default is \code{0.02}. Values below or above the default decrease
#' or increase all symbol dimensions proportionally. Relative adjustments can be 
#' specified directly, for example, \code{scale = 0.02 * 0.3} produces marker 
#' dimensions equal to 30 percent of the default.
#' @param size Character string controlling display of study-level point estimates:
#'  \describe{
#'    \item{"fisher"}{Study-symbol dimensions represent percentage contributions 
#'      to the estimated logit sensitivity and logit specificity of the study’s own subgroup.
#'      Weights sum to 100 separately within each subgroup and outcome dimension. 
#'      If \code{sensspec_constrain} is used in \code{\link{fitReitsmaSubgroupTMB}} 
#'      to constrain logit
#'      sensitivity and/or logit specificity to a common value across
#'      subgroups, weights for the constrained outcome dimension(s) instead
#'      sum to 100 across all studies and subgroups combined, reflecting
#'      that all studies then inform a single shared parameter. Default.}
#'    \item{"equal"}{All studies shown with equal size.}
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
#'       larger dimension measure across studies equals eight before adding
#'       the fixed baseline dimension.}
#'    \item{"se_revman"}{
#'       RevMan-style inverse-standard-error scaling. Horizontal and
#'       vertical symbol dimensions are based on approximate inverse standard
#'       errors of observed logit specificity and logit sensitivity,
#'       respectively. A common analysis-specific scaling constant is chosen
#'       so that the average larger dimension measure across studies equals
#'       eight before adding the fixed baseline dimension.}
#'  } 
#' @param HSROC if \code{TRUE}, the HSROC curve is added to the plot.
#'   Default is \code{FALSE}.
#' @param specrange A numeric vector of length 2 giving the range of
#'   specificities over which the HSROC curve is plotted.
#'   Defaults to \code{c(0.7, 0.995)}.
#' @param col Vector of colours used for subgroup-specific HSROC curves,
#'   study-level symbols, summary points, confidence and prediction region. 
#'   If \code{NULL}, colours are generated automatically.
#' @param nudge_legend Numeric horizontal offset for the subgroup legend.
#'   More negative values move the legend further right, outside the plotting area.
#'   Values closer to zero move it closer to the panel. Default is \code{-0.4}.
#' @param main Character string giving the main title of the plot.
#'   Defaults to \code{"Diagnostic Test Accuracy Meta-Analysis"}.
#' 
#' @param conflevel Confidence level for the confidence region. Default is \code{0.95}.
#' @param predlevel Confidence level for the prediction region. Default is \code{0.95}.
#' @param connectstudies Whether the point estimates of two subgroups 
#'   within the same study should be connected. Defaults to \code{FALSE}.
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
#'   \item Optional HSROC curve (if \code{HSROC = TRUE})
#' }
#'
#' Confidence and prediction regions are derived using the delta method
#' based on the estimated variance-covariance structure of the model.
#' 
#' @references
#' Harbord, R. M., Deeks, J. J., Egger, M., Whiting, P., & Sterne, J. A. C. (2007).
#' A unification of models for meta-analysis of diagnostic accuracy studies.
#' \emph{Biostatistics}, 8(2), 239--251.
#' \doi{10.1093/biostatistics/kxl004}
#' 
#' @return
#' No return value. Called for its side effect of producing a plot.
#'
#' @seealso \code{\link{fitReitsmaSubgroupTMB}}
#' @method plot ReitsmaSubgroupTMB
#' @importFrom grDevices adjustcolor rainbow
#' @export
plot.ReitsmaSubgroupTMB <- function(x,
                                    symbol=NULL,
                                    scale=0.02, 
                                    size=c("fisher","equal","sampsize","se","sampsize_revman","se_revman"), 
                                    main="Diagnostic Test Accuracy Meta-Analysis",
                                    col=NULL,
                                    nudge_legend=-0.4,
                                    HSROC=FALSE,
                                    specrange=c(0.7,0.995),
                                    conflevel=0.95,
                                    predlevel=0.95,
                                    connectstudies=FALSE,
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
    if(length(unique(x$data$subgroup)) != 2) {
      warning("'connectstudies=TRUE' is only recommended for two-subgroup comparisons." )
    }
  }
  size   <- match.arg(size)
  sub    <- levels(x$data$subgroup)
  nsub   <- length(sub)
  warn_unestimable_sroc_points(x$data)
  aes  <- validateSubgroupAesthetics(col=col,symbol=symbol,nsub=nsub)
  col  <- aes$col
  col2 <- aes$col2
  symbols2 <- aes$symbols2
  symb <- vector(mode="integer",length=nsub)
  # Calculations for percentage weights
  if(size=="fisher"){
    nstudy<- nrow(x$data)
    Y_pw  <- reshapeX_REIT(X=x$data)
    if(is.null(x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=2*nsub)
      col_names <- unlist(lapply(sub,function(sg) {c(paste0("mu_A.", sg),paste0("mu_B.", sg))}))
      colnames(X_pw) <- col_names
      se_index <- integer(nstudy)
      sp_index <- integer(nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup[i])
        g <- match(sg, sub)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        se_index[i] <- 2 * g - 1
        sp_index[i] <- 2 * g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("sens" %in% x$sensspec_constrain & !("spec" %in% x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=nsub+1)
      colnames(X_pw) <- c("mu_A.common",paste0("mu_B.", sub))
      se_index <- rep(1,nstudy)
      sp_index <- integer(nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup[i])
        g  <- match(sg, sub)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        sp_index[i] <- 1 + g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("spec" %in% x$sensspec_constrain & !("sens" %in% x$sensspec_constrain)){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=nsub+1)
      colnames(X_pw) <- c(paste0("mu_A.", sub),"mu_B.common")
      se_index <- integer(nstudy)
      sp_index <- rep(nsub+1,nstudy)
      for (i in seq_len(nstudy)) {
        sg <- as.character(x$data$subgroup[i])
        g  <- match(sg, sub)
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        se_index[i] <- g
        X_pw[row_se, se_index[i]] <- 1
        X_pw[row_sp, sp_index[i]] <- 1
      }
    }
    if("sens" %in% x$sensspec_constrain & "spec" %in% x$sensspec_constrain){
      X_pw <- matrix(0,nrow=2*nstudy,ncol=2)
      colnames(X_pw) <- c("mu_A.common","mu_B.common")
      se_index <- rep(1,nstudy)
      sp_index <- rep(2,nstudy)
      for (i in seq_len(nstudy)) {
        row_se <- 2 * i - 1
        row_sp <- 2 * i
        X_pw[row_se, 1] <- 1
        X_pw[row_sp, 2] <- 1
      } 
    }  
    XT_pw <- t(X_pw)
    Z <- diag(2*nstudy)
    invn <- 1/Y_pw$n
    invn[!is.finite(invn)] <- 100000000000 
    A <- diag(invn)
    ###
    rn <- rownames(x$sdreport2)
    mu_A  <- x$sdreport2[startsWith(rn,"mu_A."),"Estimate"]
    mu_B  <- x$sdreport2[startsWith(rn,"mu_B."),"Estimate"]
    names(mu_A) <- sub("^mu_A\\.","",names(mu_A))
    names(mu_B) <- sub("^mu_B\\.","",names(mu_B))
    mu_A2 <- mu_A[as.character(x$data$subgroup)]
    mu_B2 <- mu_B[as.character(x$data$subgroup)]
    p_pwTMB  <- x$sdreport$par.random
    p_pwTMB2 <- data.frame(sens=stats::plogis(mu_A2+p_pwTMB[names(p_pwTMB)=="sensu"]),
                           spec=stats::plogis(mu_B2+p_pwTMB[names(p_pwTMB)=="specu"]))
    p_pw <- as.vector(rbind(p_pwTMB2$sens,p_pwTMB2$spec))
    ###
    var_pw <- p_pw*(1-p_pw)
    B <- diag(var_pw)
    G <- matrix(0,ncol=nstudy*2,nrow=nstudy*2)
    for(i in 1:nstudy){
      j       <- 2*i
      sg      <- x$data$subgroup[i]
      s2_A.sg <- paste0("sigma2_A.",sg)
      s2_B.sg <- paste0("sigma2_B.",sg)
      s_AB.sg <- paste0("sigma_AB.",sg)
      varA    <- x$sdreport2[s2_A.sg,"Estimate"]
      varB    <- x$sdreport2[s2_B.sg,"Estimate"]
      sAB     <- x$sdreport2[s_AB.sg,"Estimate"]
      G[(j-1):j,(j-1):j] <- matrix(c(varA,sAB,sAB,varB),2,2)
    }
    #inverse of B (required later on)
    BI <- solve(B)
    # Create variance matrix for observations
    V <- (Z %*% G %*% t(Z)) + (A %*% BI)
    # invert the variance matrix
    invV <- solve(V)
    # derive the fishers information matrix
    fish <- XT_pw %*% invV %*% X_pw
    # invert Fishers information to obtain Var Beta hat
    varb  <- solve(fish)
    pct <- data.frame(sp = numeric(nstudy),
                      se = numeric(nstudy))
    # Get weights  
    for (i in seq_len(nstudy)){
      DM <- V
      DM[(i*2)-1, (i*2)-1] <- 100000000000
      DM[(i*2)-1, (i*2)] <- 0
      DM[(i*2), (i*2)-1] <- 0
      DM[(i*2), (i*2)] <- 100000000000
      invDM <- solve(DM)
      fishD <- XT_pw %*% invDM %*% X_pw
      fishI <- fish - fishD
      weight <- varb %*% fishI %*% varb
      se_idx  <- se_index[i]
      sp_idx  <- sp_index[i]
      pct$se[i] <- 100*(weight[se_idx,se_idx]/varb[se_idx,se_idx])
      pct$sp[i] <- 100*(weight[sp_idx,sp_idx]/varb[sp_idx,sp_idx])
    }
    pct$se <- pmax(0,pct$se)
    pct$sp <- pmax(0,pct$sp)
  } else {
    pct <- getWEIGHTS(xdata=x$data,size=size)
  }
  ####
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 4, 10),
      pty="s")
  ### Plot coordinate system
  plot_SESPGRID(main=main)
  # Plot study level estimates 
  for (i in seq_along(sub)){
    pointsXY(x=1-x$data$spec[x$data$subgroup==sub[i]], 
             y=x$data$sens[x$data$subgroup==sub[i]], 
             symbol = symbols2[i], 
             scale = scale*0.5,
             cex.x = pct$sp[x$data$subgroup==sub[i]],
             cex.y = pct$se[x$data$subgroup==sub[i]],
             col=col2[i])
    symb[i] <- switch(symbols2[i], rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
  }
  # Add the ROC curve
  if(HSROC==TRUE){
    for(i in seq_along(sub)){
      roc_points2 <- getROCpoints(Lambda=x$RutterGatsonis_recovered[sub[i],"Lambda"],
                                  beta=x$RutterGatsonis_recovered[sub[i],"beta"],
                                  specrange)
      points(roc_points2, type="l", lwd=2,ann=FALSE,col=col[i])
    }
  } ###
  # Add summary point
  for (i in seq_along(sub)){
    sg      <- sub[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    mean_point <- data.frame(1-x$sensspec[mu_B.sg,"Estimate"],
                             x$sensspec[mu_A.sg,"Estimate"])
    points(mean_point, col=col[i], cex=1.5, pch=15)
  }
  # Add confidence and prediction region
  for(i in seq_along(sub)){
    sg      <- sub[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg) 
    s2_A.sg <- paste0("sigma2_A.",sg)
    s2_B.sg <- paste0("sigma2_B.",sg)
    s_AB.sg <- paste0("sigma_AB.",sg)
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
                                 nstudy=sum(x$data$subgroup == sg, na.rm = TRUE),
                                 conflevel=conflevel,
                                 predlevel=predlevel)
    lines(region$conf, lty=2, lwd=2, col=col2[i])
    lines(region$pred, lty=3, lwd=2, col=col2[i])
  }
  ## Connect studies
  if(connectstudies){
    for(st in unique(x$data$study)){
      tmp <- x$data[x$data$study == st, ]
      tmp <- tmp[order(tmp$subgroup), ]
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
  if(HSROC==TRUE){
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "HSROC curve",
                      "Summary estimate",
                      conf_lab,
                      pred_lab,
                      "Data"), 
           pch = c(NA,NA,15,NA,NA,symb[1]), 
           lty = c(NA,1,NA,2,3,NA), 
           lwd = c(NA,2,NA,2,2,NA), 
           col = c(NA,"black","black","black","black","darkgray"))}
  else{
    legend("bottomright", 
           bty ="n",
           legend = c(NA,
                      "Summary estimate",
                      conf_lab,
                      pred_lab,
                      "Data"), 
           pch = c(NA,15,NA,NA,symb[1]), 
           lty = c(NA,NA,2,3,NA), 
           lwd = c(NA,NA,2,2,NA), 
           col = c(NA,"black","black","black","darkgray"))
  }
  legend("right",
         inset = c(nudge_legend, 0),
         legend = sub,
         col    = col,
         pch = symb,
         xpd = TRUE,
         cex = 1.2,
         bty = "n")
  invisible(NULL)
}



#' @rdname as_revman.dtametaTMB
#' @export
as_revman.ReitsmaSubgroupTMB <- function(x, ...) {
  
  sub  <- levels(x$data$subgroup)
  res  <- vector("list", length(sub))
  for(i in seq_along(sub)) {
    sg  <- sub[i]
    mu_A.sg <- paste0("mu_A.",sg)
    mu_B.sg <- paste0("mu_B.",sg)
    s2_A.sg <- paste0("sigma2_A.",sg)
    s2_B.sg <- paste0("sigma2_B.",sg)
    s_AB.sg <- paste0("sigma_AB.",sg)
    r_AB.sg <- paste0("rho_AB.",sg)
    ###
    mu_se  <- x$sdreport2[mu_A.sg,"Estimate"]
    mu_sp  <- x$sdreport2[mu_B.sg,"Estimate"]
    var_se <- x$sdreport2[s2_A.sg,"Estimate"]
    var_sp <- x$sdreport2[s2_B.sg,"Estimate"]
    cov_ss <- x$sdreport2[s_AB.sg,"Estimate"]
    cor_ss <- x$sdreport2[r_AB.sg,"Estimate"]
    ###
    Lambda <- x$RutterGatsonis_recovered[sg,"Lambda"]
    Theta  <- x$RutterGatsonis_recovered[sg,"Theta"]
    beta   <- x$RutterGatsonis_recovered[sg,"beta"]
    varA   <- x$RutterGatsonis_recovered[sg,"sigma2_alpha"]
    varT   <- x$RutterGatsonis_recovered[sg,"sigma2_theta"]
    ###
    seelse <- x$sdreport2[mu_A.sg,"Std. Error"]
    seelsp <- x$sdreport2[mu_B.sg,"Std. Error"]
    coves  <- x$vcov[mu_A.sg,mu_B.sg]
    nstudy <- sum(x$data$subgroup == sg, na.rm = TRUE)
    ###
    res[[i]] <- data.frame(
      Subgroup=sg,
      Externally_Calculated_Parameters=
        c(rep("HSROC model parameters",5),
          rep("Bivariate model parameters",6),
          rep("Confidence and prediction regions",4)),
      Parameter = c(
        "Lambda",
        "Theta",
        "beta",
        "Var(accuracy)",
        "Var(threshold)",
        "E(logitSe)",
        "E(logitSp)",
        "Var(logitSe)",
        "Var(logitSp)",
        "Cov(logits)",
        "Corr(logits)",
        "SE(E(logitSe))",
        "SE(E(logitSp))",
        "Cov(Es)",
        "Studies"),
      Estimate = c(Lambda,Theta,beta,varA,varT,
                   mu_se,mu_sp,var_se,var_sp,cov_ss,cor_ss,
                   seelse,seelsp,coves,nstudy))
  }
  
  do.call(rbind, res)
}
