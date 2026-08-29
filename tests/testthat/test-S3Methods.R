test_that("common generic methods work across inheritance hierarchy", {
  
  ## Reitsma ----------------------------------------------------
  
  data("anticcp", package = "dtametaTMB")
  
  reitsma <- fitReitsma(
    data  = anticcp,
    TP    = TP,
    FP    = FP,
    FN    = FN,
    TN    = TN,
    study = study
  )
  
  ## Reitsma subgroup -------------------------------------------
  
  data("RF", package = "dtametaTMB")
  
  RF2 <- RF[RF$method %in% c("LA","ELISA","Nephelometry"), ]
  RF2$method <- factor(
    RF2$method,
    levels = c("LA","ELISA","Nephelometry")
  )
  
  reitsmaSub <- fitReitsmaSubgroup(
    data     = RF2,
    TP       = TP,
    FP       = FP,
    FN       = FN,
    TN       = TN,
    study    = study,
    subgroup = method
  )
  
  ## Reitsma LCA ------------------------------------------------
  
  data("pap", package = "dtametaTMB")
  
  reitsmaLCA <- fitReitsmaLCA(
    data  = pap,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = id
  )
  
  ## Reitsma subgroup LCA ---------------------------------------
  
  reitsmaSubLCA <- fitReitsmaSubgroupLCA(
    data  = pap,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = id,
    subgroup = type
  )
  
  ## Rutter-Gatsonis --------------------------------------------
  
  ruga <- fitRutterGatsonis(
    data  = RF,
    TP    = TP,
    FP    = FP,
    FN    = FN,
    TN    = TN,
    study = study
  )
  
  ## Rutter-Gatsonis subgroup ----------------------------------
  
  rugaSub <- fitRutterGatsonisSubgroup(
    data      = RF2,
    TP        = TP,
    FP        = FP,
    FN        = FN,
    TN        = TN,
    study     = study,
    subgroup  = method,
    constrain = "shape"
  )
  
  ## Rutter-Gatsonis meta-regression ----------------------------
  Z  <- model.matrix(~method,data=RF2)
  Z2 <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  Z_pred <- matrix(c(1,0,0,1,1,0,1,0,1),ncol=3,nrow=3,byrow=T)
  constrain <- list(shape_coef=factor(c(1, rep(NA, ncol(Z2) - 1))))
  
  rugaReg <- fitRutterGatsonisReg(data=RF2,
                                  TP=TP,
                                  FP=FP,
                                  FN=FN,
                                  TN=TN,
                                  study=study,
                                  Z=Z2,
                                  Z_pred=Z_pred,
                                  map=constrain)
 
  
  ## Rutter-Gatsonis LCA ---------------------------------------
  
  rugaLCA <- fitRutterGatsonisLCA(
    data  = pap,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = id
  )
  
  ## Rutter-Gatsonis subgroup LCA -------------------------------
  
  rugaSubLCA <- fitRutterGatsonisSubgroupLCA(
    data     = pap,
    y11      = y11,
    y10      = y10,
    y01      = y01,
    y00      = y00,
    study    = id,
    subgroup = type
  )
  
  ## Hoyer ------------------------------------------------------
  
  data("diabetes", package = "dtametaTMB")
  
  hoyer <- fitHoyer(
    data            = diabetes,
    TP              = TP,
    FP              = FP,
    FN              = FN,
    TN              = TN,
    threshold       = threshold,
    study           = study,
    smallest        = 2,
    largest         = 10,
    dist            = "loglogistic",
    testdirection   = "greater",
    eval_threshold  = c(5, 6, 7)
  )
  
  models <- list(
    reitsma,
    reitsmaSub,
    reitsmaLCA,
    reitsmaSubLCA,
    ruga,
    rugaSub,
    rugaReg,
    rugaLCA,
    rugaSubLCA,
    hoyer
  )
  
  models2 <- list(
    reitsma,
    reitsmaSub,
    reitsmaLCA,
    reitsmaSubLCA,
    ruga,
    rugaSub,
    rugaLCA,
    rugaSubLCA
  )
  
  
  ## -----------------------------------------------------------
  ## Generic dispatch tests
  ## -----------------------------------------------------------
  
  for(mod in models){
    
    expect_no_error(coef(mod))
    expect_no_error(vcov(mod))
    expect_no_error(logLik(mod))
    expect_no_error(AIC(mod))
    expect_no_error(BIC(mod))
    
    expect_true(length(coef(mod)) > 0)
    
    expect_true(
      inherits(logLik(mod), "logLik")
    )
  }
  
  for(mod in models2){
    expect_no_error(as_revman(mod))
  }
  
})