test_that("plot.HoyerAFT runs without error (greater)", {
  data("diabetes", package = "dtametaTMB")
  
  res <- restructure_data(
    data = diabetes,
    TP = TP, FP = FP, FN = FN, TN = TN,
    threshold = threshold,
    study = study,
    smallest = 2,
    largest = 10,
    testdirection = "greater"
  )
  
  init <- initHoyerAFT(res$restructured, dist = "loglogistic")
  
  fit <- fitHoyerAFT(res, init)
  
  expect_no_error(plot(fit))
})



test_that("plot.HoyerAFT returns NULL invisibly", {
  data("diabetes", package = "dtametaTMB")
  
  res <- restructure_data(
    data = diabetes,
    TP = TP, FP = FP, FN = FN, TN = TN,
    threshold = threshold,
    study = study,
    smallest = 2,
    largest = 10
  )
  
  init <- initHoyerAFT(res$restructured, dist = "loglogistic")
  
  fit <- fitHoyerAFT(res, init)
  
  out <- plot(fit)
  
  expect_null(out)
})


test_that("plot.HoyerAFT works for all distributions", {
  data("diabetes", package = "dtametaTMB")
  
  res <- restructure_data(
    data = diabetes,
    TP = TP, FP = FP, FN = FN, TN = TN,
    threshold = threshold,
    study = study,
    smallest = 2,
    largest = 10
  )
  
  dists <- c("weibull", "lognormal", "loglogistic")
  
  for (d in dists) {
    init <- initHoyerAFT(res$restructured, dist = d)
    fit  <- fitHoyerAFT(res, init)
    
    expect_no_error(plot(fit))
  }
})


test_that("plot.HoyerAFT supports all size and symbol options", {
  data("diabetes")
  
  fit <- fitHoyer(
    data = diabetes,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    threshold = threshold,
    study = study,
    smallest = 2,
    largest = 10
  )
  
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  expect_no_error(plot(fit, size = "se_revman"))
  
  expect_no_error(plot(fit, symbol = "rectangle"))
  expect_no_error(plot(fit, symbol = "ellipse"))
  expect_no_error(plot(fit, symbol = "diamond"))
  expect_no_error(plot(fit, symbol = "triangle"))
  expect_no_error(plot(fit, symbol = "plus"))
  expect_no_error(plot(fit, symbol = "cross"))
  expect_no_error(plot(fit, symbol = "star"))
})
