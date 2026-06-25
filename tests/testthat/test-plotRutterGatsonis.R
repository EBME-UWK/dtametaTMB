test_that("plot.RutterGatsonis runs without error", {
  skip_on_cran()
  
  obj <- list(
    data = data.frame(
      sens = c(0.8, 0.7),
      spec = c(0.9, 0.85)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.5, 0.1),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  expect_no_error(plot(obj))
})

test_that("plot.RutterGatsonis returns NULL invisibly", {
  skip_on_cran()
  
  obj <- list(
    data = data.frame(
      sens = c(0.8, 0.7),
      spec = c(0.9, 0.85)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.5, 0.1),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  expect_null(plot(obj))
})

test_that("plot.RutterGatsonis works with multiple studies", {
  skip_on_cran()
  
  obj <- list(
    data = data.frame(
      sens = c(0.6, 0.7, 0.8),
      spec = c(0.85, 0.8, 0.9)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.4, 0.2),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  expect_no_error(plot(obj))
})

test_that("plot.RutterGatsonis handles extreme values without crashing", {
  skip_on_cran()
  
  obj <- list(
    data = data.frame(
      sens = c(0.01, 0.99),
      spec = c(0.01, 0.99)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.3, -0.2),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  expect_no_error(plot(obj))
})

test_that("plot adjusts correctly when ROC subset is empty", {
  skip_on_cran()
  
  # Force restrictive range so roc_points2 may become empty
  obj <- list(
    data = data.frame(
      sens = c(0.5),
      spec = c(0.5)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.2, 0.1),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  expect_no_error(plot(obj))
})

test_that("plot does not permanently change par settings", {
  skip_on_cran()
  
  old_par <- par(no.readonly = TRUE)
  
  obj <- list(
    data = data.frame(
      sens = c(0.8, 0.7),
      spec = c(0.9, 0.85)
    ),
    sdreport2 = data.frame(
      Estimate = c(0.5, 0.1),
      row.names = c("Lambda", "beta")
    )
  )
  
  class(obj) <- "RutterGatsonis"
  
  plot(obj)
  
  new_par <- par(no.readonly = TRUE)
  
  expect_equal(old_par$pty, new_par$pty)
})