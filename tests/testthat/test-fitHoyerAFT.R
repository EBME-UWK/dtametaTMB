test_that("fitHoyerAFT returns model object", {
  data("diabetes")
  dat <- diabetes

  res <- restructure(dat, TP, FP, FN, TN, threshold, study, 2, 10)
  init <- getInitParms(res$restructured)

  fit <- fitHoyerAFT(res, init)

  expect_s3_class(fit, "HoyerAFT")
  expect_true("sdreport2" %in% names(fit))
})

test_that("invalid threshold is rejected", {
  data("diabetes")
  dat <- diabetes

  res <- restructure(dat, TP, FP, FN, TN, threshold, study, 2, 10)
  init <- getInitParms(res$restructured)

  expect_error(
    fitHoyerAFT(res, init, threshold = -1)
  )
})



test_that("init validation fails if missing columns", {
  data_list <- list(
    original = data.frame(threshold = 1:3,testdirection=c("greater","greater","greater")),
    restructured = data.frame(
      lowerB = 1, upperB = 2,
      events0 = 1, events1 = 1,
      ctype = 1, study = 1
    )
  )
  
  bad_init <- data.frame(beta0_init = 0)
  
  expect_error(
    fitHoyerAFT(data_list, bad_init),
    "must be a valid output from getInitParms"
  )
})

test_that("threshold validation works", {
  data_list <- list(
    original = data.frame(threshold = c(2, 3, 4),
                          testdirection=c("greater","greater","greater")),
    restructured = data.frame(
      lowerB = c(1, 2),
      upperB = c(2, 3),
      events0 = c(1, 1),
      events1 = c(1, 1),
      ctype = c(1, 1),
      study = c(1, 2)
    )
  )
  
  init <- data.frame(
    beta0_init = 0,
    lambda0_init = 1,
    beta1_init = 0,
    lambda1_init = 1,
    su0_init = 0.5,
    su1_init = 0.5,
    coru0u1_init = 0,
    distcode = 1
  )
  
  expect_error(
    fitHoyerAFT(data_list, init, threshold = "a"),
    "must be numeric"
  )
  
  expect_error(
    fitHoyerAFT(data_list, init, threshold = c(1, Inf)),
    "finite"
  )
  
  expect_error(
    fitHoyerAFT(data_list, init, threshold = c(-1, 2)),
    "must be positive"
  )
})

test_that("sensspec output structure is correct", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  
  data_list <- list(
    original = data.frame(
      threshold = c(2, 3, 4),
      testdirection=c("greater","greater","greater")
    ),
    restructured = data.frame(
      lowerB = c(1, 2, 3),
      upperB = c(2, 3, 4),
      events0 = c(10, 12, 8),
      events1 = c(30, 25, 20),
      ctype = c(1, 1, 1),
      study = c(1, 2, 3)
    )
  )
  
  init <- data.frame(
    beta0_init = 0,
    lambda0_init = 1,
    beta1_init = 0,
    lambda1_init = 1,
    su0_init = 0.5,
    su1_init = 0.5,
    coru0u1_init = 0,
    distcode = 1
  )
  
  fit <- fitHoyerAFT(data_list, init)
  
  expect_true(is.data.frame(fit$sensspec))
  expect_true(all(c(
    "threshold", "conflevel",
    "Sens", "SensCI_Lower", "SensCI_Upper",
    "Spec", "SpecCI_Lower", "SpecCI_Upper"
  ) %in% colnames(fit$sensspec)))
})

test_that("custom threshold is respected", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  
  data_list <- list(
    original = data.frame(
      threshold = c(2, 3, 4),
      testdirection=c("greater","greater","greater")
    ),
    restructured = data.frame(
      lowerB = c(1, 2, 3),
      upperB = c(2, 3, 4),
      events0 = c(10, 12, 8),
      events1 = c(30, 25, 20),
      ctype = c(1, 1, 1),
      study = c(1, 2, 3)
    )
  )
  
  init <- data.frame(
    beta0_init = 0,
    lambda0_init = 1,
    beta1_init = 0,
    lambda1_init = 1,
    su0_init = 0.5,
    su1_init = 0.5,
    coru0u1_init = 0,
    distcode = 1
  )
  
  fit <- fitHoyerAFT(data_list, init, threshold = 5)
  
  expect_equal(fit$sensspec$threshold[1], 5)
})
