#' @keywords internal
#' @noRd
rename_hsrocsubLCA_rows <- function(ss, lsub) {
  
  rn <- rownames(ss)
  if (is.null(rn)) rn <- rep("", nrow(ss))
  
  i_prv   <- which(rn == "prev_coef")
  i_acc   <- which(rn == "accuracy_coef")
  i_thr   <- which(rn == "threshold_coef")
  i_shp   <- which(rn == "shape_coef")
  
  i_acc_pred <- which(rn == "Lambda_Pred")
  i_thr_pred <- which(rn == "Theta_Pred")
  i_shp_pred <- which(rn == "beta_Pred")
  i_prv_pred <- which(rn == "prev_Pred")
  
  
  i_s2p   <- which(rn == "sigma2_prev")
  
  i_mAr   <- which(rn == "mu_A_ref")
  i_mBr   <- which(rn == "mu_B_ref")
  
  prv_names   <- c(paste0("mu_prev.", lsub[1]), paste0("nu_prev.",lsub[-1]))
  acc_names   <- c(paste0("Lambda_",  lsub[1]), paste0("xi_",     lsub[-1]))
  thr_names   <- c(paste0("Theta_",   lsub[1]), paste0("gamma_",  lsub[-1]))
  shp_names   <- c(paste0("beta_",    lsub[1]), paste0("delta_",  lsub[-1]))
  
  acc_names_pred  <- c(paste0("Lambda_", lsub))
  thr_names_pred  <- c(paste0("Theta_",  lsub))
  shp_names_pred  <- c(paste0("beta_",   lsub))
  prv_names_pred  <- c(paste0("mu_prev.",lsub))
  
  s2_prev_names   <- c(paste0("sigma2_prev.", lsub))
  
  mAr_name        <- "mu_A.ref"
  mBr_name        <- "mu_B.ref"
  
  rn[i_prv] <- prv_names
  rn[i_acc] <- acc_names
  rn[i_thr] <- thr_names
  rn[i_shp] <- shp_names
  
  rn[i_acc_pred] <- acc_names_pred
  rn[i_thr_pred] <- thr_names_pred
  rn[i_shp_pred] <- shp_names_pred
  rn[i_prv_pred] <- prv_names_pred
  
  rn[i_s2p] <- s2_prev_names
  rn[i_mAr] <- mAr_name
  rn[i_mBr] <- mBr_name

  rownames(ss) <- rn
  ss
}