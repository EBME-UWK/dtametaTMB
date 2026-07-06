#' Export Model Results for RevMan
#'
#' Converts fitted model objects into a format suitable for manual entry
#' into the Diagnostic Test Accuracy module of Review Manager (RevMan).
#'
#' @param x A fitted model object.
#' @param ... Currently not used.
#'
#' @return
#' A data frame with columns:
#' \describe{
#'   \item{Externally Calculated Parameters}{Revman model/parameter types.}
#'   \item{Parameter}{RevMan parameter name.}
#'   \item{Estimate}{Parameter estimate.}
#' }
#'
#' The returned table contains the model parameters displayed by RevMan:
#' \itemize{
#'   \item \code{E(logitSe)}
#'   \item \code{E(logitSp)}
#'   \item \code{Var(logitSe)}
#'   \item \code{Var(logitSp)}
#'   \item \code{Cov(logits)}
#'   \item \code{Corr(logits)}
#'   \item \code{Lambda}
#'   \item \code{Theta}
#'   \item \code{beta}
#'   \item \code{Var(accuracy)}
#'   \item \code{Var(threshold)}
#'   \item \code{SE(E(logitSe))}
#'   \item \code{SE(E(logitSp))}
#'   \item \code{Cov(Es)}
#'   \item \code{Studies}
#' }
#'
#' @note For Rutter and Gatsonis models \code{SE(E(logitSe))}, \code{SE(E(logitSp))}, and \code{Cov(Es)} are not computed.
#' @export
#' @name as_revman.dtametaTMB
as_revman <- function(x, ...) {
  UseMethod("as_revman")
}


#' @rdname as_revman.dtametaTMB
#' @export
as_revman.Reitsma <- function(x, ...) {
  
  ## extract these from x
  mu_se   <- x$estimates[1,"Estimate"]
  mu_sp   <- x$estimates[2,"Estimate"]
  var_se  <- x$estimates[3,"Estimate"]
  var_sp  <- x$estimates[4,"Estimate"]
  cov_ss  <- x$estimates[5,"Estimate"]
  cor_ss  <- cov_ss/sqrt(var_se*var_sp)
  
  Lambda  <- x$RutterGatsonis_recovered$Lambda
  Theta   <- x$RutterGatsonis_recovered$Theta
  beta    <- x$RutterGatsonis_recovered$beta
  varA    <- x$RutterGatsonis_recovered$sigma2_alpha
  varT    <- x$RutterGatsonis_recovered$sigma2_theta
  
  seelse  <- x$estimates[1,"Std_Error"]
  seelsp  <- x$estimates[2,"Std_Error"]
  coves   <- x$vcov[1,2]
  nstudy  <- nrow(x$data)
  
  ret <- data.frame(
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
      "Studies"
    ),
    Estimate = c(Lambda,Theta,beta,varA,varT,
                 mu_se,mu_sp,var_se,var_sp,cov_ss,cor_ss,
                 seelse,seelsp,coves,nstudy), 
      
    row.names = NULL
  )
  return(ret)
}


#' @rdname as_revman.dtametaTMB
#' @export
as_revman.RutterGatsonis <- function(x, ...) {
  ## extract these from x
  mu_se   <- x$Reitsma_recovered$mu_A.sens
  mu_sp   <- x$Reitsma_recovered$mu_B.spec
  var_se  <- x$Reitsma_recovered$sigma2_A.sens
  var_sp  <- x$Reitsma_recovered$sigma2_B.spec
  cov_ss  <- x$Reitsma_recovered$sigma_AB
  cor_ss  <- cov_ss/sqrt(var_se*var_sp)
  
  Lambda  <- x$sdreport2[1,"Estimate"]
  Theta   <- x$sdreport2[2,"Estimate"]
  beta    <- x$sdreport2[3,"Estimate"]
  varA    <- x$sdreport2[4,"Estimate"]
  varT    <- x$sdreport2[5,"Estimate"]
  
  seelse  <- NA_real_
  seelsp  <- NA_real_
  coves   <- NA_real_
  nstudy  <- nrow(x$data)
  
  ret <- data.frame(
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
      "Studies"
    ),
    Estimate = c(Lambda,Theta,beta,varA,varT,
                 mu_se,mu_sp,var_se,var_sp,cov_ss,cor_ss,
                 seelse,seelsp,coves,nstudy), 
    
    row.names = NULL
  )
  return(ret)
}


#' @rdname as_revman.dtametaTMB
#' @export
as_revman.ReitsmaSubgroup <- function(x, ...) {
  
  subs     <- x$subgroups
  llsub    <- length(subs)
  scounter <- llsub*2+1
  n_study  <- table(x$data$subgroup)
  res  <- vector("list", llsub)
  for(i in seq_along(subs)) {
    sg <- subs[i]
    j  <- 2*i
    ###
    mu_se  <- x$estimates_mu[j-1,"Estimate"]
    mu_sp  <- x$estimates_mu[j,"Estimate"]
    var_se <- x$estimates_mu[scounter,"Estimate"]
    var_sp <- x$estimates_mu[scounter+1,"Estimate"]
    cov_ss <- x$estimates_mu[scounter+2,"Estimate"]
    cor_ss <- cov_ss/sqrt(var_se*var_sp)
    ###
    Lambda <- x$RutterGatsonis_recovered[sg,]$Lambda
    Theta  <- x$RutterGatsonis_recovered[sg,]$Theta
    beta   <- x$RutterGatsonis_recovered[sg,]$beta
    varA   <- x$RutterGatsonis_recovered[sg,]$sigma2_alpha
    varT   <- x$RutterGatsonis_recovered[sg,]$sigma2_theta
    ###
    seelse <- x$estimates_mu[j-1,"Std_Error"]
    seelsp <- x$estimates_mu[j,"Std_Error"]
    coves  <- x$vcov_mu[j-1,j]
    nstudy <- n_study[i]
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


#' @rdname as_revman.dtametaTMB
#' @export
as_revman.RutterGatsonisSubgroup <- function(x, ...) {
  
  subs     <- x$subgroups
  llsub    <- length(subs)
  #scounter <- llsub*2+1
  n_study  <- table(x$data$subgroup)
  res  <- vector("list", llsub)
  for(i in seq_along(subs)) {
    sg <- subs[i]
    j  <- 2*i
    ###
    mu_se  <- x$Reitsma_recovered[sg,]$mu_A.sens
    mu_sp  <- x$Reitsma_recovered[sg,]$mu_B.spec
    var_se <- x$Reitsma_recovered[sg,]$sigma2_A.sens
    var_sp <- x$Reitsma_recovered[sg,]$sigma2_B.spec
    cov_ss <- x$Reitsma_recovered[sg,]$sigma_AB
    cor_ss <- cov_ss/sqrt(var_se*var_sp)
      ###
      lamb <- paste0("Lambda_",sg)
      thet <- paste0("Theta_",sg)
      bet  <- paste0("beta_",sg)
    ###
    Lambda <- x$sdreport2[lamb,"Estimate"]
    Theta  <- x$sdreport2[thet,"Estimate"]
    beta   <- x$sdreport2[bet,"Estimate"]
    varA   <- x$sdreport2["sigma2_alpha","Estimate"]
    varT   <- x$sdreport2["sigma2_theta","Estimate"]
    ###
    seelse <- NA_real_
    seelsp <- NA_real_
    coves  <- NA_real_
    nstudy <- n_study[i]
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
