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
      "RutterGatsonis_recovered",
      "subgroups"
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







