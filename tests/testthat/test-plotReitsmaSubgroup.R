test_that("plot.ReitsmaSubgroup runs on Anti-CCP example", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation
  )
  
  expect_silent(
    plot(fit)
  )
})

test_that("plot.ReitsmaSubgroup works for unequal variances", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    variances = "unequal"
  )
  
  expect_silent(
    plot(fit)
  )
})


test_that("plot.ReitsmaSubgroup works for unequal variances", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    variances = "unequal"
  )
  
  expect_silent(
    plot(fit)
  )
})


test_that("plot.ReitsmaSubgroup accepts custom colours", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation
  )
  
  expect_silent(
    plot(
      fit,
      col = c("red", "blue")
    )
  )
})


test_that("plot.ReitsmaSubgroup connects paired studies", {
  
  fit <- fitReitsmaSubgroup(
    data = subset(schuetz, indirect == 0),
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = test
  )
  
  expect_silent(
    plot(
      fit,
      connectstudies = TRUE
    )
  )
})

test_that("plot.ReitsmaSubgroup rejects invalid size argument", {
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation
  )
  
  expect_error(
    plot(fit, size = "banana")
  )
})

