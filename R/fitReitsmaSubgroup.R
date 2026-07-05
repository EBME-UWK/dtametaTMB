#' Fit Reitsma Subgroup Model
#'
#' Fits the Reitsma bivariate random-effects model with a single categorical covariate
#' for diagnostic test accuracy (DTA) meta-analysis using a binomial-normal 
#' likelihood via \code{glmmTMB}. The model is fitted twice, once with dummy/reference
#' cell coding, and once with central-mean parameterization.
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
#' @param subgroup_constrain Optional character vector specifying
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
#'   \code{subgroup_constrain = c("sens", "spec")}, which forces
#'   all subgroup-specific sensitivity and specificity parameters
#'   to be equal across subgroups.
#'
#'   If \code{NULL} (default), separate sensitivity and specificity
#'   parameters are estimated for each subgroup.
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#'
#' @return A list of class \code{"ReitsmaSubgroup"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities
#'   \item \code{glmmTMB_mu}: fitted model object with cell means parameterization.
#'   \item \code{estimates_mu}: parameter estimates with SE with cell means parameterization.
#'   \item \code{vcov_mu}: variance-covariance matrix with cell means parameterization.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{glmmTMB_nu}: fitted model object with dummy/reference-cell parameterization
#'   \item \code{estimates_nu}: parameter estimates with SE with dummy/reference-cell parameterization.
#'   \item \code{vcov_nu}: variance-covariance matrix with dummy/reference-cell parameterization.
#'   \item \code{RutterGatsonis_recovered}: Recovered parameters in the Rutter-Gatsonis (HSROC) parameterization.
#'   \item \code{subgroups}: The subgroup levels used in the model fit.
#'   \item \code{constrain}: Random effects parameters fixed at zero.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsmaSubgroup(
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
#' @importFrom glmmTMB glmmTMB getME fixef
#' @importFrom stats binomial complete.cases qnorm plogis qlogis vcov cor sd reshape
#' @export
fitReitsmaSubgroup <- function(data,
                       TP, FP, FN, TN,
                       study,
                       subgroup,
                       constrain=NULL,
                       subgroup_constrain=NULL,
                       conflevel=0.95) {
  
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
  lsub       <- levels(X$subgroup)
  llsub      <- length(lsub)
  lsub_safe  <- make.names(lsub)
  X$subgroup_safe <- factor(make.names(X$subgroup),levels=lsub_safe)
  
  XP <- getXP(X=X)
  
  ### Get initial values
  init <- fitReitsma(data=X,
                     TP=TP,FP=FP,FN=FN,TN=TN,study=study,
                     constrain=NULL,
                     conflevel=conflevel)$estimates
  muA_init     <- init["mu_A.sens","Estimate"]
  muB_init     <- init["mu_B.spec","Estimate"]
  sA_init      <- sqrt(init["sigma2_A.sens","Estimate"])
  sA_init      <- max(sA_init,1e-05)
  sB_init      <- sqrt(init["sigma2_B.spec","Estimate"])
  sB_init      <- max(sB_init,1e-05)
  sAB_init     <- init["sigma_AB","Estimate"]
  rAB_init     <- sAB_init/(sA_init*sB_init)
  rAB_init     <- max(min(rAB_init,0.99),-0.99)
  theta3_init  <- rAB_init/sqrt(1-rAB_init**2)
  ###
  start_list_nu<- list(beta=c(muA_init,muB_init,rep(0,2*(llsub-1))),
                       theta=c(log(sA_init),log(sB_init),theta3_init))
  start_list_mu<- list(beta = rep(c(muA_init, muB_init),llsub), 
                       theta=c(log(sA_init),log(sB_init),theta3_init))
  ### Constraints
  
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
  
  
  allowed_subgroup_constraints <- c(
    "sens",
    "spec"
  )
  
  if (!is.null(subgroup_constrain)) {
    if (!is.character(subgroup_constrain)) {
      stop("'subgroup_constrain' must be a character vector or NULL.")
    }
    invalid_constraints <- setdiff(
      subgroup_constrain,
      allowed_subgroup_constraints
    )
    if(length(invalid_constraints) > 0){
      stop(
        "Unknown subgroup constraint(s): ",
        paste(invalid_constraints, collapse=", ")
      )
    }
  }
  
  map_mu <- map_nu <- list()
  
  if(!is.null(constrain)){
    if (constrain == "sigma_AB") {
      start_list_mu$theta[3] <- 0
      start_list_nu$theta[3] <- 0
      map_mu$theta <- map_nu$theta <- factor(c(1,2,NA))
    }
    if (constrain == "sigma2_A") {
      start_list_mu$theta[1] <- log(.Machine$double.eps)
      start_list_mu$theta[3] <- 0
      start_list_nu$theta[1] <- log(.Machine$double.eps)
      start_list_nu$theta[3] <- 0
      map_mu$theta <- map_nu$theta <- factor(c(NA,1,NA))
    }
    if (constrain == "sigma2_B") {
      start_list_mu$theta[2] <- log(.Machine$double.eps)
      start_list_mu$theta[3] <- 0
      start_list_nu$theta[2] <- log(.Machine$double.eps)
      start_list_nu$theta[3] <- 0
      map_mu$theta <- map_nu$theta <- factor(c(1,NA,NA))
    }
    if (constrain == "all") {
      start_list_mu$theta[1] <- log(.Machine$double.eps)
      start_list_mu$theta[2] <- log(.Machine$double.eps)
      start_list_mu$theta[3] <- 0
      start_list_nu$theta[1] <- log(.Machine$double.eps)
      start_list_nu$theta[2] <- log(.Machine$double.eps)
      start_list_nu$theta[3] <- 0
      map_mu$theta <- map_nu$theta <- factor(c(NA,NA,NA))
    }
  }
  
  ### Subgroup constraints
  if (!is.null(subgroup_constrain)) {
    beta_map_mu <- seq_len(2 * llsub)
    if ("sens" %in% subgroup_constrain) {
      beta_map_mu[seq(1, 2 * llsub, by = 2)] <- 1
    }
    if ("spec" %in% subgroup_constrain) {
      beta_map_mu[seq(2, 2 * llsub, by = 2)] <- 2
    }
      map_mu$beta <- factor(beta_map_mu)
  }
  
  if (!is.null(subgroup_constrain)){
    beta_map_nu <- seq_len(2 * llsub)
    if ("sens" %in% subgroup_constrain) {
      idxA <- seq(3, 2 * llsub, by = 2)
      start_list_nu$beta[idxA] <- 0
      beta_map_nu[idxA] <- NA
    }
    if ("spec" %in% subgroup_constrain) {
      idxB <- seq(4, 2 * llsub, by = 2)
      start_list_nu$beta[idxB] <- 0
      beta_map_nu[idxB] <- NA
    }
      map_nu$beta <- factor(beta_map_nu)
  }

  
  
  
  ### resphaping the data
  Y    <- reshapeX_REIT(X)
  ### Fitting the Reitsma model
  muA_name <- paste0("mu_A.",lsub_safe[1])
  muB_name <- paste0("mu_B.",lsub_safe[1])
  Y[[muA_name]] <- Y$sens
  Y[[muB_name]] <- Y$spec
  ###
  nuA_names <- paste0("nu_A.", lsub_safe[-1])
  nuB_names <- paste0("nu_B.", lsub_safe[-1])
  if(llsub==1){ 
    nuA_names <- nuB_names <- "" 
  } else {
    for (j in 2:llsub) {
      Y[[nuA_names[j-1]]] <- as.numeric(Y$subgroup == lsub[j]) * Y$sens
      Y[[nuB_names[j-1]]] <- as.numeric(Y$subgroup == lsub[j]) * Y$spec
    }
  }
  base_terms <- paste(c("0", muA_name, muB_name), collapse=" + ")
  nu_terms   <- c(rbind(nuA_names, nuB_names))
  fixed_rhs  <- paste(c(base_terms, nu_terms), collapse = " + ")
  form_txt   <- paste0("cbind(true, n - true) ~ ",
                        fixed_rhs,
                        " + (0 + sens + spec | recordid)")
  form       <- as.formula(form_txt)
  
  MA_Y_nu <- glmmTMB::glmmTMB(formula=form, 
                              data=Y, family=stats::binomial(link="logit"),
                              start=start_list_nu,
                              map = if (length(map_nu) == 0) NULL else map_nu)
  if (MA_Y_nu$fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  ### SAS variance covariance matrix
  theta_nu   <- glmmTMB::getME(MA_Y_nu,"theta")
  beta_fix_nu<- glmmTMB::fixef(MA_Y_nu)$cond
  V_full_nu  <- vcov(MA_Y_nu, full = TRUE)
  V_full_nu[is.na(V_full_nu)] <- 0
  esti_V_g_nu<- get2esti_V_g(beta_fix=beta_fix_nu, 
                             theta=theta_nu, 
                             V_full=V_full_nu) 
  ####
  ####
  muA_names <- paste0("mu_A.", lsub_safe)
  muB_names <- paste0("mu_B.", lsub_safe)
  for (j in 1:llsub) {
    Y[[muA_names[j]]] <- as.numeric(Y$subgroup == lsub[j]) * Y$sens
    Y[[muB_names[j]]] <- as.numeric(Y$subgroup == lsub[j]) * Y$spec
  }
  mu_terms   <- c(rbind(muA_names, muB_names))
  fixed_rhs  <- paste(mu_terms, collapse = " + ")
  form_txt   <- paste0("cbind(true, n - true) ~ 0 + ",
                       fixed_rhs,
                       " + (0 + sens + spec | recordid)")
  form       <- as.formula(form_txt)
  ####
  ####
  MA_Y_mu <- glmmTMB::glmmTMB(formula=form,
                              data=Y, family=stats::binomial(link="logit"),
                              start=start_list_mu,
                              map = if (length(map_mu) == 0) NULL else map_mu)
  if (MA_Y_mu$fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  ### SAS variance covariance matrix
  theta_mu   <- glmmTMB::getME(MA_Y_mu,"theta")
  beta_fix_mu<- glmmTMB::fixef(MA_Y_mu)$cond
  V_full_mu  <- vcov(MA_Y_mu, full = TRUE)
  V_full_mu[is.na(V_full_mu)] <- 0
  esti_V_g_mu<- get2esti_V_g(beta_fix=beta_fix_mu, 
                             theta=theta_mu, 
                             V_full=V_full_mu) 
  ### Sensitivity and Specificity
  ma_Y_mu        <- summary(MA_Y_mu)
  qq             <- stats::qnorm(1-(1-conflevel)/2)
  sesp           <- as.data.frame(ma_Y_mu$coefficients$cond)
  sesp$Orig      <- with(sesp,stats::plogis(Estimate))
  sesp$conflevel <- conflevel
  sesp$CI_Lower  <- with(sesp,stats::plogis(Estimate-qq*`Std. Error`))
  sesp$CI_Upper  <- with(sesp,stats::plogis(Estimate+qq*`Std. Error`))
  sesp           <- sesp[,(5:8)]
  sesp$type      <- sub("^mu_([AB])\\..*$", "\\1", rownames(sesp))
  sesp$type      <- c(A = "sens", B = "spec")[sesp$type]
  sesp           <- sesp[,c("type","Orig","conflevel","CI_Lower","CI_Upper")]
  ### Recover HSROC parameters
  ruga2 <- data.frame()
  scounter <- llsub*2+1
  for(i in 1:llsub){
    j <- 2*i
    ruga <- getRUGA(lsens=esti_V_g_mu$esti[j-1,"Estimate"],
                    lspec=esti_V_g_mu$esti[j  ,"Estimate"],
                    sigma_a=sqrt(esti_V_g_mu$esti[scounter,"Estimate"]),
                    sigma_b=sqrt(esti_V_g_mu$esti[scounter+1,"Estimate"]),
                    sigma_ab=esti_V_g_mu$esti[scounter+2,"Estimate"])
    ruga2 <- rbind(ruga2,ruga)
  }
  rownames(ruga2) <- lsub
  
  ##
  ret <- list(data         = XP,
              glmmTMB_mu   = MA_Y_mu,
              estimates_mu = esti_V_g_mu$esti,
              vcov_mu      = esti_V_g_mu$V_g,
              sensspec     = sesp,
              glmmTMB_nu   = MA_Y_nu,
              estimates_nu = esti_V_g_nu$esti,
              vcov_nu      = esti_V_g_nu$V_g,
              RutterGatsonis_recovered = ruga2,
              constrain    = constrain,
              subgroups    = lsub)
  class(ret) <- c("ReitsmaSubgroup","CochraneSubgroup")
  return(ret)
}