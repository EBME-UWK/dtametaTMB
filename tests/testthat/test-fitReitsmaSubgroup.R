test_that("fitReitsmaSubgroup returns a ReitsmaSubgroup object", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_s3_class(fit, "ReitsmaSubgroup")
})

test_that("fitReitsmaSubgroup returns expected components", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_named(
    fit,
    c(
      "data",
      "glmmTMB_mu",
      "estimates_mu",
      "vcov_mu",
      "sensspec",
      "glmmTMB_nu",
      "estimates_nu",
      "vcov_nu",
      "LRDOR",
      "RutterGatsonis_recovered",
      "constrain",
      "subgroups",
      "variances"
    )
  )
})

test_that("all subgroup levels are recovered", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_equal(
    length(fit$subgroups),
    nlevels(factor(RF$method))
  )
})

test_that("sensspec contains expected columns", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_true(is.data.frame(fit$sensspec))
  
  expect_true(
    all(
      c(
        "type",
        "Orig",
        "conflevel",
        "CI_Lower",
        "CI_Upper"
      ) %in% names(fit$sensspec)
    )
  )
})

test_that("sensspec estimates lie between 0 and 1", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_true(all(fit$sensspec$Orig >= 0))
  expect_true(all(fit$sensspec$Orig <= 1))
})

test_that("HSROC parameters are recovered for every subgroup", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_equal(
    nrow(fit$RutterGatsonis_recovered),
    length(fit$subgroups)
  )
})

test_that("logLik, AIC and BIC work for ReitsmaSubgroup", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_s3_class(
    logLik(fit),
    "logLik"
  )
  
  expect_true(is.numeric(AIC(fit)))
  expect_true(is.numeric(BIC(fit)))
})


test_that("summary.ReitsmaSubgroup returns expected structure", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  s <- summary(fit)
  
  expect_type(s, "list")
  
  expect_named(
    s,
    c(
      "estimates",
      "sensspec",
      "RutterGatsonis_recovered",
      "subgroups"
    )
  )
})


test_that("fitReitsmaSubgroup rejects non-data.frame input", {
  
  expect_error(
    fitReitsmaSubgroup(
      data = matrix(1:4, 2, 2),
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = subgroup
    ),
    "'data' must be a data.frame"
  )
})

test_that("all subgroup levels are recovered", {
  
  fit <- fitReitsmaSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_equal(
    nrow(fit$RutterGatsonis_recovered),
    length(fit$subgroups)
  )
  
  expect_equal(
    length(fit$subgroups),
    nlevels(factor(RF$method))
  )
})




test_that("fitReitsmaSubgroup validates constrain argument", {
  
  expect_error(
    fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      constrain = "banana"
    ),
    "constrain"
  )
  
  expect_error(
    fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      constrain = c("sigma_AB", "sigma2_A")
    ),
    "constrain"
  )
  
})


test_that("fitReitsmaSubgroup sigma_AB fixes covariance", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    constrain = "sigma_AB"
  )
  
  expect_equal(
    fit$estimates_mu["sigma_AB", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
  expect_equal(
    fit$estimates_nu["sigma_AB", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})

test_that("fitReitsmaSubgroup sigma2_A fixes sensitivity variance", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    constrain = "sigma2_A"
  )
  
  expect_lt(
    fit$estimates_mu["sigma2_A.sens", "Estimate"],
    1e-12
  )
  
  expect_lt(
    fit$estimates_nu["sigma2_A.sens", "Estimate"],
    1e-12
  )
  
  expect_equal(
    fit$estimates_mu["sigma_AB", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})

test_that("fitReitsmaSubgroup sigma2_B fixes specificity variance", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    constrain = "sigma2_B"
  )
  
  expect_lt(
    fit$estimates_mu["sigma2_B.spec", "Estimate"],
    1e-12
  )
  
  expect_lt(
    fit$estimates_nu["sigma2_B.spec", "Estimate"],
    1e-12
  )
  
  expect_equal(
    fit$estimates_mu["sigma_AB", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})

test_that("fitReitsmaSubgroup all fixes all random-effects parameters", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    constrain = "all"
  )
  
  expect_lt(
    fit$estimates_mu["sigma2_A.sens", "Estimate"],
    1e-12
  )
  
  expect_lt(
    fit$estimates_mu["sigma2_B.spec", "Estimate"],
    1e-12
  )
  
  expect_equal(
    fit$estimates_mu["sigma_AB", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})

test_that("fitReitsmaSubgroup still estimates subgroup means", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    constrain = "all"
  )
  
  expect_true(
    all(
      is.finite(
        fit$estimates_mu$Estimate
      )
    )
  )
  
})



test_that("RutterGatsonis_recovered correctly round-trips per subgroup (catches sens/spec swaps)", {
  
  data("anticcp")
  
  fit <- fitReitsmaSubgroup(
    data     = anticcp,
    TP       = TP,
    FP       = FP,
    FN       = FN,
    TN       = TN,
    study    = study,
    subgroup = generation
  )
  
  lsub      <- fit$subgroups
  lsub_safe <- make.names(lsub)
  
  expect_equal(nrow(fit$RutterGatsonis_recovered), length(lsub))
  expect_equal(rownames(fit$RutterGatsonis_recovered), lsub)
  
  for (i in seq_along(lsub)) {
    
    # HSROC-space estimates reported for this subgroup
    hsroc_row <- fit$RutterGatsonis_recovered[i, ]
    
    # Independently invert them back to Reitsma-space using the
    # package's own (separately-tested) inverse transform
    recovered <- getREIT(
      Lambda       = hsroc_row$Lambda,
      Theta        = hsroc_row$Theta,
      beta         = hsroc_row$beta,
      sigma2_alpha = hsroc_row$sigma2_alpha,
      sigma2_theta = hsroc_row$sigma2_theta
    )
    
    # Original Reitsma-space estimates for this subgroup, as
    # actually fitted by fitReitsmaSubgroup()
    muA_name <- paste0("mu_A.", lsub_safe[i])
    muB_name <- paste0("mu_B.", lsub_safe[i])
    
    original_muA <- fit$estimates_mu[muA_name, "Estimate"]
    original_muB <- fit$estimates_mu[muB_name, "Estimate"]
    original_sigma2_A <- fit$estimates_mu["sigma2_A.sens", "Estimate"]
    original_sigma2_B <- fit$estimates_mu["sigma2_B.spec", "Estimate"]
    original_sigma_AB <- fit$estimates_mu["sigma_AB", "Estimate"]
    
    # Round-trip must recover the original subgroup-specific estimates
    expect_equal(recovered$mu_A.sens, original_muA,
                 tolerance = 1e-6,
                 label = paste("mu_A round-trip for subgroup", lsub[i]))
    expect_equal(recovered$mu_B.spec, original_muB,
                 tolerance = 1e-6,
                 label = paste("mu_B round-trip for subgroup", lsub[i]))
    expect_equal(recovered$sigma2_A.sens, original_sigma2_A,
                 tolerance = 1e-6,
                 label = paste("sigma2_A round-trip for subgroup", lsub[i]))
    expect_equal(recovered$sigma2_B.spec, original_sigma2_B,
                 tolerance = 1e-6,
                 label = paste("sigma2_B round-trip for subgroup", lsub[i]))
    expect_equal(recovered$sigma_AB, original_sigma_AB,
                 tolerance = 1e-6,
                 label = paste("sigma_AB round-trip for subgroup", lsub[i]))
  }
})

