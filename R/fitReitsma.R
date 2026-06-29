#' Fit Reitsma Model
#'
#' Fits the Reitsma bivariate random-effects model for diagnostic test accuracy (DTA)
#' meta-analysis using a binomial-normal likelihood via \code{glmmTMB}.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#'
#' @return A list of class \code{"Reitsma"} with components:
#' \itemize{
#'   \item \code{data}: the original data set with derived quantities
#'   \item \code{glmmTMB}: fitted model object.
#'   \item \code{estimates}: parameter estimates with SE.
#'   \item \code{vcov}: variance-covariance matrix.
#'   \item \code{sensspec}: sensitivity and specificity estimates.
#'   \item \code{LRDOR}: DOR and likelihood ratios.
#'   \item \code{RutterGatsonis_recovered}: HSROC parameters.
#' }
#'
#' @examples
#' data("anticcp")
#' fit <- fitReitsma(
#'   data = anticcp,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study
#' )
#' fit$estimates
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
fitReitsma <- function(data,
                       TP, FP, FN, TN,
                       study, conflevel=0.95) {
  
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
      FN = data[[FN_col]]
    )
    
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
  sA_init      <- stats::sd(logit_sens)
  sA_init      <- max(sA_init,1e-05)
  sB_init      <- stats::sd(logit_spec)
  sB_init      <- max(sB_init,1e-05)
  rAB_init     <- max(min(cor(logit_sens,logit_spec),0.99),-0.99)
  theta3_init  <- rAB_init/sqrt(1-rAB_init**2)
  
  ### Reshaping the data
  X$true1 <- X$TP
  X$true0 <- X$TN 
  X$n1    <- X$TP+X$FN
  X$n0    <- X$FP+X$TN
  X$recordid <- 1:nrow(X)
  Y <- reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1", "true0")), 
               timevar="sens", times=c(1,0), v.names=c("n","true")) 
  Y <- Y[order(Y$id),]  
  Y$spec <- 1-Y$sens
  ### Fitting the Reitsma model
  MA_Y <- glmmTMB::glmmTMB(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec | id), 
                           data=Y, family=stats::binomial(link="logit"),
                           start=list(beta=c(muA_init,muB_init),
                                      theta=c(log(sA_init),log(sB_init),theta3_init)))
  if (MA_Y$fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  ma_Y <- summary(MA_Y)
  S         <- ma_Y$vcov$cond
  qq        <- stats::qnorm(1-(1-conflevel)/2)
  res       <- as.data.frame(ma_Y$coefficients$cond)
  res$Orig  <- with(res,plogis(Estimate))
  res$conflevel <- conflevel
  res$CI_Lower  <- with(res,plogis(Estimate-qq*`Std. Error`))
  res$CI_Upper  <- with(res,plogis(Estimate+qq*`Std. Error`))
  res       <- res[,(5:8)];
  colnames(res) <- c("Estimate","conflevel","CI_Lower","CI_Upper")
  ### SAS variance covariance matrix
  theta     = glmmTMB::getME(MA_Y,"theta")
  beta_fix  = glmmTMB::fixef(MA_Y)$cond
  V_full    = vcov(MA_Y, full = TRUE)
  ### Define the transformed parameter vector
  g <- c(mu_A = beta_fix[1],
         mu_B = beta_fix[2],
         sigma2_A.sens = exp(2*theta[1]),
         sigma2_B.spec = exp(2*theta[2]),
         sigma_AB = (theta[3] / sqrt(1 + theta[3]**2)) * exp(theta[1]) * exp(theta[2]))
  ### Implement Jacobian
  J = matrix(0, nrow = 5, ncol = 5)
  ### fixed effects
  J[1, 1]  = 1
  J[2, 2]  = 1
  ### variances
  J[3, 3]  = 2 * exp(2 * theta[1])
  J[4, 4]  = 2 * exp(2 * theta[2])
  ### covariance
  sigma_A  = exp(theta[1])
  sigma_B  = exp(theta[2])
  rho      = theta[3] / sqrt(1 + theta[3]**2)
  sigma_AB = rho * sigma_A * sigma_B
  J[5, 3]  = sigma_AB # d/d theta1
  J[5, 4]  = sigma_AB # d/d theta2
  J[5, 5]  = sigma_A * sigma_B / (1 + theta[3]**2)**(3/2)
  rownames(J) = colnames(J) = c("mu_A.sens",
                                "mu_B.spec",
                                "sigma2_A.sens",
                                "sigma2_B.spec",
                                "sigma_AB")
  ### Apply delta method to get SAS variance-covariance matrix
  V_g  = J %*% V_full %*% t(J)
  rep <- data.frame("Estimate"=g,"Std_Error"=sqrt(diag(V_g)))
  #rep$CI_Lower <- with(rep,Estimate-qq*Std_Error)
  #rep$CI_Upper <- with(rep,Estimate+qq*Std_Error)

  # diagnostic odds ratio, the positive and negative
  # likelihood ratios
  lsens  <- rep[1,"Estimate"]
  lspec  <- rep[2,"Estimate"]
  DOR    <- exp(lsens+lspec) 
  LRp    <- plogis(lsens)/(1-plogis(lspec))
  LRn    <- ((1-plogis(lsens))/plogis(lspec)) 

  se.logDOR = as.numeric(sqrt(c(1,1) %*% S %*% c(1,1)))
  dLRp = rbind(1/(1+exp(lsens)),exp(lspec)/(1+exp(lspec)))
  se.logLRp = as.numeric(sqrt(t(dLRp) %*% S %*% dLRp))
  dLRn = rbind(exp(lsens)/(1+exp(lsens)),1/(1+exp(lspec)))
  se.logLRn = as.numeric(sqrt(t(dLRn) %*% S %*% dLRn))

  rep2 <- data.frame(Estimate = c(DOR, LRp, LRn), 
                     conflevel=conflevel,
                     CI_Lower = c(exp(log(DOR)-qq*se.logDOR), exp(log(LRp)-qq*se.logLRp), exp(log(LRn)-qq*se.logLRn)), 
                     CI_Upper = c(exp(log(DOR)+qq*se.logDOR), exp(log(LRp)+qq*se.logLRp), exp(log(LRn)+qq*se.logLRn)),
                     row.names = c("DOR", "LR+", "LR-")) 
  # Recover Rutter and Gatsonis estimates
  sigma2_a <- rep[3,"Estimate"]
  sigma2_b <- rep[4,"Estimate"]
  sigma_ab <- rep[5,"Estimate"]
  sigma_a  <- sqrt(sigma2_a)
  sigma_b  <- sqrt(sigma2_b)
  rep3 <- data.frame(
            Lambda   = (((sigma_b/sigma_a)**0.5) * lsens) + ((sigma_a/sigma_b)**0.5 *lspec),
            Theta    = 0.5*((((sigma_b/sigma_a)**0.5 )*lsens) - (((sigma_a/sigma_b)**0.5) *lspec)),
            beta     = log(sigma_b/sigma_a),
            sigma2_alpha = 2*((sigma_a*sigma_b) + sigma_ab),
            sigma2_theta = 0.5*((sigma_a*sigma_b) - sigma_ab),
            row.names="Estimate (recovered)")
  ret <- list(data=XP,
              glmmTMB=MA_Y,
              estimates=rep,
              vcov=V_g,
              sensspec=res,
              LRDOR=rep2,
              RutterGatsonis_recovered=rep3)
  class(ret) <- c("Reitsma","Cochrane")
  return(ret)
}