test_that("plot.RutterGatsonisSubgroup runs without error", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(fit)
  )
})

test_that("plot.RutterGatsonisSubgroup accepts user colours", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      col = c("black", "red", "blue", "green")
    )
  )
})


test_that("all size options work", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
})



test_that("custom specrange works", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      specrange = c(0.8, 0.99)
    )
  )
})


test_that("custom title works", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      main = "Test Plot"
    )
  )
})

test_that("plot dispatches correctly", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_true(
    is.function(
    getS3method(
      "plot",
      "RutterGatsonisSubgroup"
    )
    )
  )
})

test_that("plot.RutterGatsonisSubgroup runs without error", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(fit)
  )
})

test_that("plot.RutterGatsonisSubgroup supports all size options", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
})

test_that("plot.RutterGatsonisSubgroup accepts custom colours", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      col = c("red", "blue", "green", "black")
    )
  )
})

test_that("plot.RutterGatsonisSubgroup accepts custom specrange", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      specrange = c(0.8, 0.99)
    )
  )
})


test_that("plot.RutterGatsonisSubgroup accepts graphical options", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_no_error(
    plot(
      fit,
      main = "My Plot",
      nudge_legend = -0.2
    )
  )
})

test_that("plot.RutterGatsonisSubgroup rejects invalid size", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  expect_error(
    plot(
      fit,
      size = "banana"
    )
  )
})

test_that("plot dispatches correctly for RutterGatsonisSubgroup", {
  
  expect_true(
    is.function(
    getS3method(
      "plot",
      "RutterGatsonisSubgroup"
    )
    )
  )
})

test_that("plot.RutterGatsonisSubgroup restores graphics parameters", {
  
  fit <- fitRutterGatsonisSubgroup(
    data = RF,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = method
  )
  
  oldpar <- par(no.readonly = TRUE)
  
  expect_no_error(plot(fit))
  
  expect_equal(
    par("mar"),
    oldpar$mar
  )
})






