test_that("common generic methods work with edge case data set", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  
  edge <- data.frame(study=1:5,
                     subgroup=c("A","A","A","B","B"),
                     TP=0,FP=0,FN=0,TN=0,
                     y11=0,y10=0,y01=0,y00=0)
  edge[1,"y11"] <- edge[1,"TP"] <- 1
  edge[4,"y10"] <- edge[4,"FP"] <- 2
  edge[1,"y10"] <- edge[1,"FP"] <- 1
  edge[4,"y01"] <- edge[4,"FN"] <- 2
  
  # edge <- tub
  # colnames(edge)[1] <- "study"
  # edge$subgroup <- c(rep(c("A","B"),14),"A")
  # edge$TP <- edge$y11
  # edge$FP <- edge$y10
  # edge$FN <- edge$y01
  # edge$TN <- edge$y00
  
  ## Reitsma ----------------------------------------------------
  
  reitsma <- fitReitsma(
    data  = edge,
    TP    = TP,
    FP    = FP,
    FN    = FN,
    TN    = TN,
    study = study
  )
  
  reitsmaTMB <- fitReitsmaTMB(
    data  = edge,
    TP    = TP,
    FP    = FP,
    FN    = FN,
    TN    = TN,
    study = study
  )
  
  
  ## Reitsma subgroup -------------------------------------------
  
  reitsmaSub <- fitReitsmaSubgroup(
    data     = edge,
    TP       = TP,
    FP       = FP,
    FN       = FN,
    TN       = TN,
    study    = study,
    subgroup = subgroup
  )
  
  reitsmaSubTMB <- fitReitsmaSubgroupTMB(
    data     = edge,
    TP       = TP,
    FP       = FP,
    FN       = FN,
    TN       = TN,
    study    = study,
    subgroup = subgroup
  )
  
  ## Reitsma LCA ------------------------------------------------
  
  
  reitsmaLCA <- fitReitsmaLCA(
    data  = edge,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = study
  )
  
  ## Reitsma subgroup LCA ---------------------------------------
  
  reitsmaSubLCA <- fitReitsmaSubgroupLCA(
    data  = edge,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = study,
    subgroup = subgroup
  )
  
  ## Rutter-Gatsonis --------------------------------------------
  
  ruga <- fitRutterGatsonis(
    data  = edge,
    TP    = TP,
    FP    = FP,
    FN    = FN,
    TN    = TN,
    study = study
  )
  
  ## Rutter-Gatsonis subgroup ----------------------------------
  
  rugaSub <- fitRutterGatsonisSubgroup(
    data      = edge,
    TP        = TP,
    FP        = FP,
    FN        = FN,
    TN        = TN,
    study     = study,
    subgroup  = subgroup
  )
  
  ## Rutter-Gatsonis LCA ---------------------------------------
  
  rugaLCA <- fitRutterGatsonisLCA(
    data  = edge,
    y11   = y11,
    y10   = y10,
    y01   = y01,
    y00   = y00,
    study = study
  )
  
  ## Rutter-Gatsonis subgroup LCA -------------------------------
  
  rugaSubLCA <- fitRutterGatsonisSubgroupLCA(
    data     = edge,
    y11      = y11,
    y10      = y10,
    y01      = y01,
    y00      = y00,
    study    = study,
    subgroup = subgroup
  )
  
  ## Hoyer ------------------------------------------------------
  edge[4,"TP"] <- edge[4,"TN"] <- edge[4,"y11"] <- edge[4,"y00"] <- 4
  
  hoyer <- fitHoyer(
    data            = edge,
    TP              = TP,
    FP              = FP,
    FN              = FN,
    TN              = TN,
    threshold       = study,
    study           = study,
    smallest        = 0.5,
    largest         = 30,
    dist            = "loglogistic",
    testdirection   = "greater")
  
  models <- list(
    reitsma,
    reitsmaTMB,
    reitsmaSub,
    reitsmaSubTMB,
    reitsmaLCA,
    reitsmaSubLCA,
    ruga,
    rugaSub,
    rugaLCA,
    rugaSubLCA,
    hoyer
  )
  
  models2 <- list(
    reitsma,
    reitsmaTMB,
    reitsmaSub,
    reitsmaSubTMB,
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
    expect_no_error(print(mod))
    expect_no_error(summary(mod))
    expect_no_error(print(summary(mod)))
    expect_no_error(plot(mod))
    expect_no_error(plot(mod,studyCI=TRUE,studylabels=TRUE))
    expect_no_error(forest(mod))
  }
  
  for(mod in models2){
    expect_no_error(as_revman(mod))
  }
  
  for(mod in models2){
    expect_no_error(plot(mod,studyCI=TRUE,HSROC=TRUE))
  }

})