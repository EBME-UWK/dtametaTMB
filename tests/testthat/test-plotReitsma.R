test_that("plot.Reitsma runs without error", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:5,
    TP = c(30, 25, 40, 35, 20),
    FP = c(10, 15, 5, 8, 12),
    FN = c(5, 10, 8, 7, 6),
    TN = c(50, 45, 60, 55, 48)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_no_error(plot(fit))
})

test_that("plot.Reitsma returns NULL invisibly", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_null(plot(fit))
})

test_that("plot.Reitsma works with HSROC = TRUE", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:5,
    TP = c(30, 25, 40, 35, 20),
    FP = c(10, 15, 5, 8, 12),
    FN = c(5, 10, 8, 7, 6),
    TN = c(50, 45, 60, 55, 48)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_no_error(plot(fit, HSROC = TRUE))
})

test_that("plot.Reitsma respects scale parameter", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:5,
    TP = c(30, 25, 40, 35, 20),
    FP = c(10, 15, 5, 8, 12),
    FN = c(5, 10, 8, 7, 6),
    TN = c(50, 45, 60, 55, 48)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_no_error(plot(fit, scale = 0.005))
  expect_no_error(plot(fit, scale = 0.02))
})

test_that("plot.Reitsma handles extreme values", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:4,
    TP = c(1, 99, 40, 35),
    FP = c(99, 1, 5, 8),
    FN = c(1, 99, 8, 7),
    TN = c(99, 1, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  expect_no_error(plot(fit))
})

test_that("plot.Reitsma does not permanently change par settings", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  
  dat <- data.frame(
    study = 1:4,
    TP = c(30, 25, 40, 35),
    FP = c(10, 15, 5, 8),
    FN = c(5, 10, 8, 7),
    TN = c(50, 45, 60, 55)
  )
  
  fit <- fitReitsma(dat, TP, FP, FN, TN, study)
  
  old_par <- par(no.readonly = TRUE)
  plot(fit)
  new_par <- par(no.readonly = TRUE)
  
  expect_equal(old_par$pty, new_par$pty)
})


test_that("plot.Reitsma errors on malformed input", {
  obj <- list()
  class(obj) <- "Reitsma"
  
  expect_error(plot(obj))
})