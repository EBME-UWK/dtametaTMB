test_that("plot.RutterGatsonisSubgroup runs without error", {
  skip_on_cran()
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
  skip_on_cran()
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
      col = c("black", "red", "blue", "green","orange")
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
  
  expect_no_error(plot(fit, size = "fisher"))
  expect_no_error(plot(fit, size = "fisher_revman"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  expect_no_error(plot(fit, size = "se_revman"))
})



test_that("custom specrange works", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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


test_that("plot.RutterGatsonisSubgroup accepts graphical options", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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



