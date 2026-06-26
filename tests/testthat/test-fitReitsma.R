test_that("fitReitsma runs with data frame input", {
  dat <- data.frame(
    study = 1:5,
    TP = c(30, 25, 40, 35, 20),
    FP = c(10, 15, 5, 8, 12),
    FN = c(5, 10, 8, 7, 6),
    TN = c(50, 45, 60, 55, 48)
  )
  
  fit <- fitReitsma(
    data = dat,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study
  )
  
  expect_s3_class(fit, "Reitsma")
  expect_type(fit, "list")
  
  expect_named(fit, c(
    "data", "glmmTMB", "estimates", "vcov",
    "sensspec", "LRDOR", "RutterGatsonis_recovered"
  ))
})

test_that("non-numeric input throws error", {
  dat <- data.frame(
    study = 1:3,
    TP = c("x", "y", "z"),
    FP = c(1, 2, 3),
    FN = c(1, 2, 3),
    TN = c(1, 2, 3)
  )
  
  expect_error(
    fitReitsma(dat, TP, FP, FN, TN, study),
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
    fitReitsma(dat_neg, TP, FP, FN, TN, study),
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
    fitReitsma(dat_nonint, TP, FP, FN, TN, study),
    "must contain non-negative integer counts"
  )
})

test_that("estimates output structure is correct", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_true(is.data.frame(fit$estimates))
  expect_true(all(c("Estimate", "Std_Error") %in% colnames(fit$estimates)))
})

test_that("sensspec output structure is correct", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_true(is.data.frame(fit$sensspec))
  expect_true(all(c("Estimate", "conflevel", "CI_Lower", "CI_Upper") %in% colnames(fit$sensspec)))
})

test_that("LRDOR output contains DOR and likelihood ratios", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_true(is.data.frame(fit$LRDOR))
  expect_true(all(c("DOR", "LR+", "LR-") %in% rownames(fit$LRDOR)))
})

test_that("variance-covariance matrix is symmetric", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_true(is.matrix(fit$vcov))
  expect_equal(fit$vcov, t(fit$vcov), tolerance = 1e-8)
})

test_that("RutterGatsonis recovered parameters exist", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  expect_true(is.data.frame(fit$RutterGatsonis_recovered))
  expect_true(all(c(
    "Lambda", "Theta", "beta",
    "sigma2_alpha", "sigma2_theta"
  ) %in% colnames(fit$RutterGatsonis_recovered)))
})

test_that("alpha affects CI width in estimates", {
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit1 <- fitReitsma(dat, TP, FP, FN, TN, study, conflevel = 0.95)
  fit2 <- fitReitsma(dat, TP, FP, FN, TN, study, conflevel = 0.80)
  
  width1 <- fit1$estimates$CI_Upper - fit1$estimates$CI_Lower
  width2 <- fit2$estimates$CI_Upper - fit2$estimates$CI_Lower
  
  expect_true(all(width2 < width1))
})