#' @keywords internal
#' @noRd
rename_reitsubLCA_rows <- function(ss, lsub) {
  
  rn <- rownames(ss)
  if (is.null(rn)) rn <- rep("", nrow(ss))
  
  i_nup   <- which(rn == "nu_prev")
  i_nuA   <- which(rn == "nu_A_index")
  i_nuB   <- which(rn == "nu_B_index")
  i_mAr   <- which(rn == "mu_A_ref")
  i_mBr   <- which(rn == "mu_B_ref")
  
  nu_prev_names  <- paste0("nu_prev.", lsub[-1])
  nu_A_names     <- paste0("nu_A.index.", lsub[-1])
  nu_B_names     <- paste0("nu_B.index.", lsub[-1])
  mAr_name       <- "mu_A.ref"
  mBr_name       <- "mu_B.ref"
  
  i_mup   <- which(rn == "mu_prev")
  i_muA   <- which(rn == "mu_A_index")
  i_muB   <- which(rn == "mu_B_index")
  i_s2p   <- which(rn == "sigma2_prev")
  i_s2A   <- which(rn == "sigma2_A_index")
  i_s2B   <- which(rn == "sigma2_B_index")
  i_sAB   <- which(rn == "sigma_AB_index")
  i_rAB   <- which(rn == "rho_AB_index")
  
  mu_prev_names     <- paste0("mu_prev.", lsub)
  mu_A_names        <- paste0("mu_A.index.", lsub)
  mu_B_names        <- paste0("mu_B.index.", lsub)
  sigma2_prev_names <- paste0("sigma2_prev.", lsub)
  sigma2_A_names    <- paste0("sigma2_A.index.", lsub)
  sigma2_B_names    <- paste0("sigma2_B.index.", lsub)
  sigma_AB_names    <- paste0("sigma_AB.index.", lsub)
  rho_AB_names      <- paste0("rho_AB.index.", lsub)
  
  rn[i_nup] <- nu_prev_names
  rn[i_nuA] <- nu_A_names
  rn[i_nuB] <- nu_B_names
  
  rn[i_mup] <- mu_prev_names
  rn[i_muA] <- mu_A_names
  rn[i_muB] <- mu_B_names
  rn[i_s2p] <- sigma2_prev_names
  rn[i_s2A] <- sigma2_A_names
  rn[i_s2B] <- sigma2_B_names
  rn[i_sAB] <- sigma_AB_names
  rn[i_rAB] <- rho_AB_names
  
  rn[i_mAr] <- mAr_name
  rn[i_mBr] <- mBr_name
  
  rownames(ss) <- rn
  ss
}