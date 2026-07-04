test_that("summary.RutterGatsonisSubgroup returns expected structure", {
  
  fit <- fitRutterGatsonisSubgroup(
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
      "Reitsma_recovered",
      "subgroups"
    )
  )
})


test_that("summary.RutterGatsonisSubgroup returns object components correctly", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  s <- summary(fit)
  
  expect_equal(s$estimates, fit$sdreport2)
  expect_equal(s$sensspec, fit$sensspec)
  expect_equal(s$Reitsma_recovered, fit$Reitsma_recovered)
  expect_equal(s$subgroups, fit$subgroups)
})


test_that("summary.RutterGatsonisSubgroup components have expected classes", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  s <- summary(fit)
  
  
  expect_true(is.data.frame(s$sensspec))
  expect_true(is.data.frame(s$Reitsma_recovered))
  
  expect_true(
    is.character(s$subgroups) ||
      is.factor(s$subgroups)
  )
})

test_that("summary.RutterGatsonisSubgroup contains subgroup results", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  s <- summary(fit)
  
  expect_gt(nrow(s$estimates), 0)
  expect_gt(nrow(s$sensspec), 0)
  expect_gt(length(s$subgroups), 1)
})


test_that("fitRutterGatsonisSubgroup validates constraints", {
  
  expect_error(
    fitRutterGatsonisSubgroup(
      data = RF,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = method,
      constrain = "banana"
    ),
    "Unknown constraint"
  )
  
  expect_error(
    fitRutterGatsonisSubgroup(
      data = RF,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = method,
      constrain = c("shape", "banana")
    ),
    "Unknown constraint"
  )
  
})

test_that("sigma2_alpha fixes accuracy heterogeneity", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "sigma2_alpha"
  )
  
  expect_lt(
    fit$sdreport$value["sigma2_alpha"],
    1e-12
  )
  
})


test_that("sigma2_theta fixes threshold heterogeneity", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "sigma2_theta"
  )
  
  expect_lt(
    fit$sdreport$value["sigma2_theta"],
    1e-12
  )
  
})


test_that("shape_zero fixes all shape parameters", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "shape_zero"
  )
  
  beta_rows <- grep(
    "^beta_",
    rownames(fit$sdreport2)
  )
  
  expect_true(
    all(
      abs(
        fit$sdreport2[beta_rows, "Estimate"]
      ) < 1e-12
    )
  )
  
})


test_that("shape imposes common shape parameter", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "shape"
  )
  
  beta_rows <- grep(
    "^beta_",
    rownames(fit$sdreport2)
  )
  
  expect_equal(
    length(
      unique(
        round(
          fit$sdreport2[beta_rows, "Estimate"],
          8
        )
      )
    ),
    1
  )
  
})


test_that("accuracy constraint imposes common Lambda", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "accuracy"
  )
  
  lambda_rows <- grep(
    "^Lambda_",
    rownames(fit$sdreport2)
  )
  
  expect_equal(
    length(
      unique(
        round(
          fit$sdreport2[lambda_rows, "Estimate"],
          8
        )
      )
    ),
    1
  )
  
})

test_that("threshold constraint imposes common Theta", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = "threshold"
  )
  
  theta_rows <- grep(
    "^Theta_",
    rownames(fit$sdreport2)
  )
  
  expect_equal(
    length(
      unique(
        round(
          fit$sdreport2[theta_rows, "Estimate"],
          8
        )
      )
    ),
    1
  )
  
})

test_that("multiple subgroup constraints work together", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method,
    constrain = c(
      "accuracy",
      "shape",
      "sigma2_alpha"
    )
  )
  
  expect_lt(
    fit$sdreport$value["sigma2_alpha"],
    1e-12
  )
  
})


