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

