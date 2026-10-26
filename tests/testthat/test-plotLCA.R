## ReitsmaLCA

test_that("default ReitsmaLCA plot works", {
  
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot returns invisible NULL", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_null(
    plot(fit)
  )
  
})

test_that("all size options work", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(plot(fit, size = "eb"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "eb_revman"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  
})


test_that("HSROC overlay works", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    plot(fit, HSROC = TRUE)
  )
  
})

test_that("custom specificity range works", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    
    plot(
      fit,
      HSROC = TRUE,
      specrange = c(0.5, 0.99)
    )
    
  )
  
})


test_that("custom confidence and prediction levels work", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    
    plot(
      fit,
      conflevel = 0.90,
      predlevel = 0.90
    )
    
  )
  
})

test_that("conflevel validation works", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_error(
    plot(fit, conflevel = 0),
    "conflevel"
  )
  
  expect_error(
    plot(fit, conflevel = 1),
    "conflevel"
  )
  
})


test_that("predlevel validation works", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_error(
    plot(fit, predlevel = 0),
    "predlevel"
  )
  
  expect_error(
    plot(fit, predlevel = 1),
    "predlevel"
  )
  
})


test_that("plot works for all supported constraints", {
  skip_on_cran()
  constraints <- list(
    NULL,
    "sigma_AB.index",
    "sigma2_A.index",
    "sigma2_B.index",
    "all"
  )
  
  fits <- lapply(
    constraints,
    function(con) {
      
      fitReitsmaLCA(
        data = pap,
        y11 = y11,
        y10 = y10,
        y01 = y01,
        y00 = y00,
        study = id,
        constrain = con
      )
      
    }
  )
  
  expect_true(
    all(
      vapply(
        fits,
        function(fit) {
          !inherits(
            try(plot(fit), silent = TRUE),
            "try-error"
          )
        },
        logical(1)
      )
    )
  )
  
})

test_that("plot does not modify fitted object", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  fit_before <- fit
  
  plot(fit)
  
  expect_identical(
    fit,
    fit_before
  )
  
})

test_that("plot works for fixed effects model", {
  skip_on_cran()
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "all"
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

## RutterGatsonisLCA

test_that("default RutterGatsonisLCA plot works", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot returns invisible NULL", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_null(
    plot(fit)
  )
  
})

test_that("all size options work", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(plot(fit, size = "eb"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "eb_revman"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  
})

test_that("custom specificity range works", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    
    plot(
      fit,
      specrange = c(0.5, 0.99)
    )
    
  )
  
})

test_that("custom title works", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_no_error(
    
    plot(
      fit,
      main = "Custom HSROC Plot"
    )
    
  )
  
})

test_that("plot works for shape constrained model", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "shape"
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot works for sigma2_alpha constrained model", {
  skip_on_cran()
  fit <- suppressWarnings(
    fitRutterGatsonisLCA(
      data = pap,
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id,
      constrain = "sigma2_alpha"
    )
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot does not modify fitted object", {
  skip_on_cran()
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  fit_before <- fit
  
  plot(fit)
  
  expect_identical(
    fit,
    fit_before
  )
  
})







## ReitsmaSubgroupLCA

test_that("default Reitsma subgroup plot works", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot returns invisible NULL", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_null(
    plot(fit)
  )
  
})

test_that("all size options work", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_no_error(plot(fit,size="eb"))
  expect_no_error(plot(fit,size="equal"))
  expect_no_error(plot(fit,size="sampsize"))
  expect_no_error(plot(fit,size="eb_revman"))
  expect_no_error(plot(fit,size="sampsize_revman"))
  
})

test_that("HSROC overlay works", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_no_error(
    plot(fit, HSROC = TRUE)
  )
  
})

test_that("custom colours work", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_no_error(
    
    plot(
      fit,
      col = c("red","blue")
    )
    
  )
  
})

test_that("conflevel validation works", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_error(
    plot(fit, conflevel = 0),
    "conflevel"
  )
  
  expect_error(
    plot(fit, conflevel = 1),
    "conflevel"
  )
  
})


test_that("predlevel validation works", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_error(
    plot(fit, predlevel = 0),
    "predlevel"
  )
  
  expect_error(
    plot(fit, predlevel = 1),
    "predlevel"
  )
  
})


test_that("plot works for fixed effects model", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    constrain = "all"
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot works with unequal variances", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    variances = "unequal"
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("plot does not modify fitted object", {
  skip_on_cran()
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  fit_before <- fit
  
  plot(fit)
  
  expect_identical(
    fit,
    fit_before
  )
  
})


## RuterGatsonisSubgroupLCA

test_that("default subgroup HSROC plot runs", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("all size options work", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(plot(fit, size = "eb"))
  expect_no_error(plot(fit, size = "equal"))
  expect_no_error(plot(fit, size = "sampsize"))
  expect_no_error(plot(fit, size = "eb_revman"))
  expect_no_error(plot(fit, size = "sampsize_revman"))
  
})


test_that("custom colours work", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(
    
    plot(
      fit,
      col = c("red", "blue")
    )
    
  )
  
})


test_that("custom specificity range works", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(
    
    plot(
      fit,
      specrange = c(0.5, 0.99)
    )
    
  )
  
})

test_that("connectstudies works for two subgroups", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(
    
    plot(
      fit,
      connectstudies = TRUE
    )
    
  )
  
})

test_that("connectstudies warns with more than two subgroups", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method
  )
  
  expect_warning(
    
    plot(
      fit,
      connectstudies = TRUE
    ),
    
    "only recommended for two-subgroup comparisons"
    
  )
  
})


test_that("plot returns invisible NULL", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_null(
    plot(fit)
  )
  
})

test_that("plot does not modify fitted object", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  fit_before <- fit
  
  plot(fit)
  
  expect_identical(
    fit,
    fit_before
  )
  
})

test_that("plot works for shape_zero constrained models", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method,
    constrain = "shape_zero"
  )
  
  expect_no_error(
    plot(fit)
  )
  
})

test_that("custom title works", {
  skip_on_cran()
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_no_error(
    
    plot(
      fit,
      main = "My HSROC Plot"
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
  fit1 <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "sens"
  )
  
  fit2 <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "spec"
  )
  
  fit3 <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = c("sens","spec")
  )
  
  expect_no_error(plot(fit1, size = "eb_revman"))
  expect_no_error(plot(fit2, size = "eb_revman"))
  expect_no_error(plot(fit3, size = "eb_revman"))
})