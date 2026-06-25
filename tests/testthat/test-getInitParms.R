test_that("getInitParms returns valid output", {
  data("diabetes")
  dat <- diabetes
  res <- restructure(dat, TP, FP, FN, TN, threshold, study, 2, 10)

  init <- getInitParms(res$restructured)

  expect_true(is.data.frame(init))
  expect_true("distcode" %in% names(init))
})


test_that("invalid distribution is rejected", {
  data("diabetes")
  dat <- diabetes

  res <- restructure(dat, TP, FP, FN, TN, threshold, study, 2, 10)

  expect_error(
    getInitParms(res$restructured, dist = "wrong")
  )
})
