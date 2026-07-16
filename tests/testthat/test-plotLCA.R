## ReitsmaLCA

test_that("default ReitsmaLCA plot works", {
  
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
  
})


test_that("HSROC overlay works", {
  
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
  
})

test_that("custom specificity range works", {
  
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
  
})

test_that("HSROC overlay works", {
  
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
  
})


test_that("custom colours work", {
  
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