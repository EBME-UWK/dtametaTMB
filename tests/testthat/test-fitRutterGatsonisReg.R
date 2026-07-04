test_that("fitRutterGatsonisReg returns a RutterGatsonisReg object", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z
  )
  
  expect_s3_class(fit, "RutterGatsonisReg")
})

test_that("fitRutterGatsonisReg returns expected components", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z
  )
  
  expect_named(
    fit,
    c(
      "data",
      "fit",
      "sdreport",
      "sdreport2",
      "sensspec"
    )
  )
})



test_that("summary.RutterGatsonisReg returns expected components", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z
  )
  
  s <- summary(fit)
  
  expect_type(s, "list")
  expect_true(length(s) > 0)
})

test_that("RutterGatsonisReg produces sensitivity estimates", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z
  )
  
  expect_true(is.data.frame(fit$sensspec))
  
  expect_true(
    all(
      c(
        "spec",
        "Sens",
        "SensCI_Lower",
        "SensCI_Upper"
      ) %in% names(fit$sensspec)
    )
  )
})


test_that("RutterGatsonisReg fits with shape = FALSE", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z)
  
  expect_equal(fit$fit$convergence, 0)
})

test_that("RutterGatsonisReg fits with shape = TRUE", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z,
  )
  
  expect_equal(fit$fit$convergence, 0)
})


test_that("fitRutterGatsonisReg rejects invalid Z rows", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  
  expect_error(
    fitRutterGatsonisReg(
      data = RF,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      Z = Z
    ),
    "two rows per study"
  )
})


test_that("fitRutterGatsonisReg rejects incompatible Z_pred", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  Z_pred <- matrix(0, nrow = 1, ncol = ncol(Z) + 1)
  
  expect_error(
    fitRutterGatsonisReg(
      data = RF,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      Z = Z,
      Z_pred = Z_pred
    ),
    "same number of columns"
  )
})


test_that("logLik, AIC and BIC work for RutterGatsonisReg", {
  
  data("RF")
  
  Z <- model.matrix(~ method, data = RF)
  Z <- Z[rep(seq_len(nrow(Z)), each = 2), , drop = FALSE]
  
  fit <- fitRutterGatsonisReg(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    Z = Z
  )
  
  expect_s3_class(logLik(fit), "logLik")
  
  expect_true(is.numeric(AIC(fit)))
  expect_true(is.numeric(BIC(fit)))
})








