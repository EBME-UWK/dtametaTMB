valid_base <- data.frame(
  study = c(1,1,1),
  threshold = c(0.1, 0.2, 0.3),
  TP = c(60,50,40),
  FP = c(20,10,5),
  FN = c(10,20,30),
  TN = c(80,90,95)
)

test_that("restructure works with default testdirection (greater)", {
  
  res <- restructure_data(
    valid_base,
    TP=TP, FP=FP, FN=FN, TN=TN,
    threshold=threshold,
    study=study,
    smallest=0.01,
    largest=0.5
  )
  
  expect_true(is.list(res))
  expect_true(is.data.frame(res$restructured))
  expect_equal(unique(res$original$testdirection), "greater")
})

test_that("restructure works for testdirection = 'less'", {
  
  dat_less <- data.frame(
    study = c(1,1,1),
    threshold = c(0.1, 0.2, 0.3),
    TP = c(40,50,60),   # increasing!
    FP = c(5,10,20),    # increasing!
    FN = c(30,20,10),   # decreasing
    TN = c(95,90,80)    # decreasing
  )
  
  res <- restructure_data(
    dat_less,
    TP=TP, FP=FP, FN=FN, TN=TN,
    threshold=threshold,
    study=study,
    smallest=0.01,
    largest=0.5,
    testdirection="less"
  )
  
  expect_true(is.list(res))
  expect_true(is.data.frame(res$restructured))
  expect_equal(unique(res$original$testdirection), "less")
})


test_that("less direction correctly swaps TP/FN and FP/TN internally", {
  
  dat_less <- data.frame(
    study = c(1,1,1),
    threshold = c(0.1,0.2,0.3),
    TP = c(40,50,60),
    FP = c(5,10,20),
    FN = c(30,20,10),
    TN = c(95,90,80)
  )
  
  res <- restructure_data(
    dat_less,
    TP=TP, FP=FP, FN=FN, TN=TN,
    threshold=threshold,
    study=study,
    smallest=0.01,
    largest=0.5,
    testdirection="less"
  )
  
  # original unchanged
  expect_equal(res$original$TP, dat_less$TP)
  
  # check totals unchanged
  expect_equal(res$original$TP + res$original$FN,
               dat_less$TP + dat_less$FN)
})

test_that("data valid for less fails for greater", {
  
  dat_less <- data.frame(
    study = c(1,1,1),
    threshold = c(0.1,0.2,0.3),
    TP = c(40,50,60),
    FP = c(5,10,20),
    FN = c(30,20,10),
    TN = c(95,90,80)
  )
  
  expect_error(
    restructure_data(
      dat_less,
      TP=TP, FP=FP, FN=FN, TN=TN,
      threshold=threshold,
      study=study,
      smallest=0.01,
      largest=0.5,
      testdirection="greater"
    )
  )
})

test_that("greater and less preserve totals", {
  
  res1 <- restructure_data(
    valid_base,
    TP=TP, FP=FP, FN=FN, TN=TN,
    threshold=threshold,
    study=study,
    smallest=0.01,
    largest=0.5,
    testdirection="greater"
  )
  
  # construct reversed version
  dat_less <- valid_base
  tmp_TP <- dat_less$TP
  tmp_FP <- dat_less$FP
  
  dat_less$TP <- dat_less$FN
  dat_less$FP <- dat_less$TN
  dat_less$FN <- tmp_TP
  dat_less$TN <- tmp_FP
  
  res2 <- restructure_data(
    dat_less,
    TP=TP, FP=FP, FN=FN, TN=TN,
    threshold=threshold,
    study=study,
    smallest=0.01,
    largest=0.5,
    testdirection="less"
  )
  
  expect_equal(res1$original$D, res2$original$D)
  expect_equal(res1$original$H, res2$original$H)
})

