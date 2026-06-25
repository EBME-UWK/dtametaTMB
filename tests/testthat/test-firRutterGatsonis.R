test_that("fitRutterGatsonis runs with data frame input", {
  dat <- data.frame(
    study = 1:5,
    TP = c(30, 25, 40, 35, 20),
    FP = c(10, 15, 5, 8, 12),
    FN = c(5, 10, 8, 7, 6),
    TN = c(50, 45, 60, 55, 48)
  )
  
  fit <- fitRutterGatsonis(
    data = dat,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study
  )
  
  expect_s3_class(fit, "RutterGatsonis")
  expect_type(fit, "list")
  expect_named(fit, c(
    "data", "fit", "sdreport", "sdreport2",
    "sensspec", "Reitsma_recovered"
  ))
})

test_that("non-numeric input throws error", {
  dat <- data.frame(
    study = 1:3,
    TP = c("a", "b", "c"),
    FP = c(1, 2, 3),
    FN = c(1, 2, 3),
    TN = c(1, 2, 3)
  )
  
  expect_error(
    fitRutterGatsonis(
      data = dat,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study
    ),
    "Columns must be numeric"
  )
})

test_that("negative or non-integer counts throw error", {
  dat_neg <- data.frame(
    study = 1:3,
    TP = c(10, -1, 5),
    FP = c(1, 2, 3),
    FN = c(1, 2, 3),
    TN = c(1, 2, 3)
  )
  
  expect_error(
    fitRutterGatsonis(dat_neg, TP, FP, FN, TN, study),
    "must contain non-negative integer counts"
  )
  
  dat_nonint <- data.frame(
    study = 1:3,
    TP = c(10.5, 2, 3),
    FP = c(1, 2, 3),
    FN = c(1, 2, 3),
    TN = c(1, 2, 3)
  )
  
  expect_error(
    fitRutterGatsonis(dat_nonint, TP, FP, FN, TN, study),
    "must contain non-negative integer counts"
  )
})

test_that("sensspec output structure is correct", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitRutterGatsonis(dat, TP, FP, FN, TN, study)
  
  expect_true(is.data.frame(fit$sensspec))
  expect_true(all(c(
    "spec", "conflevel", "logitsens",
    "Std_Error", "CI_Lower", "CI_Upper",
    "Sens", "SensCI_Lower"
  ) %in% names(fit$sensspec)))
})

test_that("Reitsma recovered parameters exist", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitRutterGatsonis(dat, TP, FP, FN, TN, study)
  
  expect_true(is.data.frame(fit$Reitsma_recovered))
  expect_true(all(c("mu_A.sens", 
                    "mu_B.spec",
                    "sigma2_A.sens",
                    "sigma2_B.spec",
                    "sigma_AB") %in% colnames(fit$Reitsma_recovered)))
})

test_that("alpha parameter affects CI width", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit1 <- fitRutterGatsonis(dat, TP, FP, FN, TN, study, conflevel = 0.95)
  fit2 <- fitRutterGatsonis(dat, TP, FP, FN, TN, study, conflevel = 0.80)
  
  width1 <- fit1$sensspec$CI_Upper - fit1$sensspec$CI_Lower
  width2 <- fit2$sensspec$CI_Upper - fit2$sensspec$CI_Lower
  
  expect_true(all(width2 < width1))  
})

test_that("custom spec is respected", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitRutterGatsonis(dat, TP, FP, FN, TN, study, spec = 0.8)
  
  expect_equal(fit$sensspec$spec[1], 0.8)
})