test_that("plot.ReitsmaSubgroup runs on Anti-CCP example", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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

test_that("all size options work", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation
  )

  expect_no_error(plot(fit, size = "fisher"))
  expect_no_error(plot(fit, size = "fisher_revman"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  expect_no_error(plot(fit, size = "se_revman"))
})


test_that("plot.ReitsmaSubgroup connects paired studies", {
  skip_on_cran()
  
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
      symbol=c("diamond","triangle"),
      size="sampsize_revman",
      scale=0.02*1.5,
      connectstudies = TRUE
    )
  )
})

test_that("plot.ReitsmaSubgroup rejects invalid size argument", {
  skip_on_cran()
  
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


test_that("all sensspec_constrain options work", {
  skip_on_cran()
  
  fit1 <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "sens"
  )
  
  fit2 <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "spec"
  )
  
  fit3 <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = c("sens","spec")
  )
  
  expect_no_error(plot(fit1, size = "fisher"))
  expect_no_error(plot(fit2, size = "fisher"))
  expect_no_error(plot(fit3, size = "fisher"))
})

###
###
test_that("plot.ReitsmaSubgroupTMB runs on Anti-CCP example", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
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


test_that("plot.ReitsmaSubgroupTMB works for unequal variances", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
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


test_that("plot.ReitsmaSubgroupTMB accepts custom colours", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
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

test_that("all size options work", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation
  )
  
  expect_no_error(plot(fit, size = "fisher"))
  expect_no_error(plot(fit, size = "fisher_revman"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "se"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  expect_no_error(plot(fit, size = "se_revman"))
})


test_that("plot.ReitsmaSubgroupTMB connects paired studies", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
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
      symbol=c("diamond","triangle"),
      scale=0.04,
      connectstudies = TRUE
    )
  )
})

test_that("plot.ReitsmaSubgroupTMB rejects invalid size argument", {
  skip_on_cran()
  
  fit <- fitReitsmaSubgroupTMB(
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


test_that("all sensspec_constrain options work", {
  skip_on_cran()
  
  fit1 <- fitReitsmaSubgroupTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "sens"
  )
  
  fit2 <- fitReitsmaSubgroupTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "spec"
  )
  
  fit3 <- fitReitsmaSubgroupTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = c("sens","spec")
  )
  
  expect_no_error(plot(fit1, size = "fisher_revman"))
  expect_no_error(plot(fit2, size = "fisher_revman"))
  expect_no_error(plot(fit3, size = "fisher_revman"))
})
