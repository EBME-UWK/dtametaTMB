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

