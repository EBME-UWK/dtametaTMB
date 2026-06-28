test_that("initHoyer returns valid output", {
  data("diabetes")
  dat <- diabetes
  res <- restructure_data(dat, TP, FP, FN, TN, threshold, study, 2, 10)

  init <- initHoyerAFT(res$restructured)

  expect_true(is.data.frame(init))
  expect_true("distcode" %in% names(init))
})


test_that("initHoyer distribution is rejected", {
  data("diabetes")
  dat <- diabetes

  res <- restructure_data(dat, TP, FP, FN, TN, threshold, study, 2, 10)

  expect_error(
    initHoyerAFT(res$restructured, dist = "wrong")
  )
})
