test_that("fitReitsmaLCA rejects invalid constrain arguments", {
  
  expect_error(
    fitReitsmaLCA(
      data = pap,
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id,
      constrain = "nonsense"
    ),
    "must be one of"
  )
  
})

test_that("fitReitsmaLCA requires a data.frame", {
  
  expect_error(
    fitReitsmaLCA(
      data = matrix(1:10, ncol = 2),
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id
    ),
    "data.frame"
  )
  
})

test_that("fitReitsmaLCA returns expected object structure", {
  
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_s3_class(fit, "ReitsmaLCA")
  expect_true("ReitsmaLCA" %in% class(fit))
  
  expect_true(is.list(fit))
  
  expect_true("sdreport2" %in% names(fit))
  expect_true("sensspec" %in% names(fit))
  expect_true("LRDOR" %in% names(fit))
  expect_true("vcov" %in% names(fit))
  
})

test_that("pap example reproduces benchmark estimates", {
  
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  est <- fit$sdreport2
  
  expect_equal(
    est["mu_prev", "Estimate"],
    0.556,
    tolerance = 0.02
  )
  
  expect_equal(
    est["mu_A.index", "Estimate"],
    0.640,
    tolerance = 0.02
  )
  
  expect_equal(
    est["mu_B.index", "Estimate"],
    1.620,
    tolerance = 0.02
  )
  
})

test_that("all sensitivity and specificity estimates lie within bounds", {
  
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  ss <- fit$sensspec
  
  expect_true(all(ss$Estimate >= 0))
  expect_true(all(ss$Estimate <= 1))
  
  expect_true(all(ss$CI_Lower >= 0))
  expect_true(all(ss$CI_Upper <= 1))
  
})


test_that("all supported constraints fit successfully", {
  
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
  
  expect_length(fits, 5)
  
  expect_true(
    all(vapply(fits, inherits, logical(1), "ReitsmaLCA"))
  )
  
})

test_that("sigma_AB.index constraint fixes covariance", {
  
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "sigma_AB.index"
  )
  
  expect_equal(
    fit$sdreport2["sigma_AB.index", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})


test_that("all constraint removes heterogeneity parameters", {
  
  fit <- fitReitsmaLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "all"
  )
  
  expect_equal(
    fit$sdreport2["sigma2_A.index", "Estimate"],
    0,
    tolerance = 1e-6
  )
  
  expect_equal(
    fit$sdreport2["sigma2_B.index", "Estimate"],
    0,
    tolerance = 1e-6
  )
  
})

test_that("rows with missing values are removed", {
  
  dat <- pap
  
  dat$y11[1] <- NA
  
  expect_message(
    
    fit <- fitReitsmaLCA(
      data = dat,
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id
    ),
    
    "Removed rows with missing values"
  )
  
  expect_equal(
    nrow(fit$data),
    nrow(pap) - 1
  )
  
})


test_that("fitRutterGatsonisLCA rejects invalid constraints", {
  
  expect_error(
    fitRutterGatsonisLCA(
      data = pap,
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id,
      constrain = "foobar"
    ),
    "Unknown constraint"
  )
  
})

test_that("fitRutterGatsonisLCA requires a data.frame", {
  
  expect_error(
    fitRutterGatsonisLCA(
      data = matrix(1:10, nrow = 5),
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id
    ),
    "data.frame"
  )
  
})

test_that("fitRutterGatsonisLCA returns expected object", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_s3_class(fit, "RutterGatsonisLCA")
  expect_true("RutterGatsonisLCA" %in% class(fit))
  
  expect_true("sdreport2" %in% names(fit))
  expect_true("sensspec" %in% names(fit))
  expect_true("Reitsma_recovered" %in% names(fit))
  
})

test_that("recovered Reitsma parameters are finite", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  rec <- fit$Reitsma_recovered
  
  expect_true(
    all(is.finite(unlist(rec)))
  )
  
})


test_that("estimated sensitivity is within bounds", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  ss <- fit$sensspec
  
  expect_true(all(ss$Sens >= 0))
  expect_true(all(ss$Sens <= 1))
  
  expect_true(all(ss$SensCI_Lower >= 0))
  expect_true(all(ss$SensCI_Upper <= 1))
  
})

test_that("default specificity equals 0.8", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id
  )
  
  expect_equal(
    fit$sensspec$spec[1],
    0.8
  )
  
})

test_that("user supplied specificity is respected", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    spec = 0.9
  )
  
  expect_equal(
    fit$sensspec$spec[1],
    0.9
  )
  
})

test_that("all supported constraints fit", {
  
  constraints <- list(
    NULL,
    "sigma2_alpha",
    "sigma2_theta",
    "shape",
    c("sigma2_alpha", "shape")
  )
  
  fits <- lapply(
    constraints,
    function(con) {
      fitRutterGatsonisLCA(
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
  
  expect_length(fits, 5)
  
  expect_true(
    all(
      vapply(
        fits,
        inherits,
        logical(1),
        "RutterGatsonisLCA"
      )
    )
  )
  
})


test_that("shape constraint fixes beta to zero", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "shape"
  )
  
  expect_equal(
    fit$sdreport2["beta", "Estimate"],
    0,
    tolerance = 1e-8
  )
  
})


test_that("sigma2_alpha constraint fixes variance", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "sigma2_alpha"
  )
  
  expect_equal(
    fit$sdreport2["sigma2_alpha", "Estimate"],
    0,
    tolerance = 1e-6
  )
  
})

test_that("sigma2_theta constraint fixes variance", {
  
  fit <- fitRutterGatsonisLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    constrain = "sigma2_theta"
  )
  
  expect_equal(
    fit$sdreport2["sigma2_theta", "Estimate"],
    0,
    tolerance = 1e-6
  )
  
})

test_that("missing rows are removed", {
  
  dat <- pap
  dat$y11[1] <- NA
  
  expect_message(
    
    fit <- fitRutterGatsonisLCA(
      data = dat,
      y11 = y11,
      y10 = y10,
      y01 = y01,
      y00 = y00,
      study = id
    ),
    
    "Removed rows with missing values"
  )
  
  expect_equal(
    nrow(fit$data),
    nrow(pap) - 1
  )
  
})



######################
### Subgroup tests ###
######################

test_that("fitReitsmaSubgroupLCA returns expected object", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_s3_class(fit, "ReitsmaSubgroupLCA")
  expect_true("ReitsmaSubgroupLCA" %in% class(fit))
  
  expect_true("subgroups" %in% names(fit))
  expect_true("sensspec" %in% names(fit))
  expect_true("LRDOR" %in% names(fit))
  expect_true("RutterGatsonis_recovered" %in% names(fit))
  
})


test_that("subgroup labels are preserved", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation
  )
  
  expect_equal(
    sort(fit$subgroups),
    sort(levels(factor(anticcp$generation)))
  )
  
})


test_that("sens constraint enforces equal sensitivities", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "sens"
  )
  
  est <- fit$sdreport2
  
  sens <- est[
    grep("^mu_A.index", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    max(sens),
    min(sens),
    tolerance = 1e-8
  )
  
})

test_that("spec constraint enforces equal specificities", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = "spec"
  )
  
  est <- fit$sdreport2
  
  spec <- est[
    grep("^mu_B.index", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    max(spec),
    min(spec),
    tolerance = 1e-8
  )
  
})

test_that("sens and spec constraints collapse subgroup means", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = anticcp,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = generation,
    sensspec_constrain = c("sens","spec")
  )
  
  est <- fit$sdreport2
  
  sens <- est[
    grep("^mu_A.index", rownames(est)),
    "Estimate"
  ]
  
  spec <- est[
    grep("^mu_B.index", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(max(sens), min(sens), tolerance = 1e-8)
  expect_equal(max(spec), min(spec), tolerance = 1e-8)
  
})


test_that("unequal variance model fits", {
  
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
  
  expect_s3_class(
    fit,
    "ReitsmaSubgroupLCA"
  )
  
})

test_that("constraints not allowed with unequal variances", {
  
  expect_error(
    
    fitReitsmaSubgroupLCA(
      data = anticcp,
      y11 = TP,
      y10 = FP,
      y01 = FN,
      y00 = TN,
      study = study,
      subgroup = generation,
      variances = "unequal",
      constrain = "all"
    ),
    
    "currently not supported"
  )
  
})

test_that("all common variance constraints fit", {
  
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
      
      fitReitsmaSubgroupLCA(
        data = anticcp,
        y11 = TP,
        y10 = FP,
        y01 = FN,
        y00 = TN,
        study = study,
        subgroup = generation,
        constrain = con
      )
      
    }
  )
  
  expect_length(fits, 5)
  
})

test_that("pap subgroup benchmark reproduced", {
  
  fit <- fitReitsmaSubgroupLCA(
    data = pap,
    y11 = y11,
    y10 = y10,
    y01 = y01,
    y00 = y00,
    study = id,
    subgroup = type
  )
  
  expect_equal(
    length(fit$subgroups),
    2
  )
  
  expect_true(
    all(fit$sensspec$Estimate > 0)
  )
  
})



test_that("fitRutterGatsonisSubgroupLCA returns expected object", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method
  )
  
  expect_s3_class(
    fit,
    "RutterGatsonisSubgroupLCA"
  )
  
  expect_true(
    "RutterGatsonisSubgroupLCA" %in% class(fit)
  )
  
  expect_true("sensspec" %in% names(fit))
  expect_true("Reitsma_recovered" %in% names(fit))
  expect_true("subgroups" %in% names(fit))
  
})



test_that("subgroup levels are preserved", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method
  )
  
  expect_equal(
    sort(fit$subgroups),
    sort(levels(factor(RF$method)))
  )
  
})


test_that("shape constraint enforces common beta", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method,
    constrain = "shape"
  )
  
  est <- fit$sdreport2
  
  beta <- est[
    grep("^beta_", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    max(beta),
    min(beta),
    tolerance = 1e-8
  )
  
})


test_that("accuracy constraint removes subgroup effects", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method,
    constrain = "accuracy"
  )
  
  est <- fit$sdreport2
  
  Lambda <- est[
    grep("^Lambda_", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    max(Lambda),
    min(Lambda),
    tolerance = 1e-8
  )
  
})

test_that("threshold constraint removes subgroup effects", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = RF,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = method,
    constrain = "threshold"
  )
  
  est <- fit$sdreport2
  
  Theta <- est[
    grep("^Theta_", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    max(Theta),
    min(Theta),
    tolerance = 1e-8
  )
  
})


test_that("sigma2_alpha constrained model fits", {
  
  expect_no_error(
    
    fitRutterGatsonisSubgroupLCA(
      data = RF,
      y11 = TP,
      y10 = FP,
      y01 = FN,
      y00 = TN,
      study = study,
      subgroup = method,
      constrain = "sigma2_alpha"
    )
    
  )
  
})

test_that("sigma2_theta constrained model fits", {
  
  expect_no_error(
    
    fitRutterGatsonisSubgroupLCA(
      data = RF,
      y11 = TP,
      y10 = FP,
      y01 = FN,
      y00 = TN,
      study = study,
      subgroup = method,
      constrain = "sigma2_theta"
    )
    
  )
  
})

test_that("shape_zero fixes all beta parameters", {
  
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
  
  est <- fit$sdreport2
  
  beta <- est[
    grep("^beta_", rownames(est)),
    "Estimate"
  ]
  
  expect_equal(
    unname(beta),
    rep(0, length(beta)),
    tolerance = 1e-8
  )
  
  
})

test_that("all major constraint combinations fit", {
  
  constraints <- list(
    NULL,
    "shape",
    "shape_zero",
    "accuracy",
    "threshold",
    c("shape","accuracy"),
    c("shape","threshold"),
    c("sigma2_alpha","shape")
  )
  
  fits <- lapply(
    constraints,
    function(con){
      
      fitRutterGatsonisSubgroupLCA(
        data = RF,
        y11 = TP,
        y10 = FP,
        y01 = FN,
        y00 = TN,
        study = study,
        subgroup = method,
        constrain = con
      )
      
    }
  )
  
  expect_length(
    fits,
    length(constraints)
  )
  
})


test_that("schuetz produces two subgroup estimates", {
  
  fit <- fitRutterGatsonisSubgroupLCA(
    data = schuetz,
    y11 = TP,
    y10 = FP,
    y01 = FN,
    y00 = TN,
    study = study,
    subgroup = test
  )
  
  expect_equal(
    length(fit$subgroups),
    2
  )
  
  expect_equal(
    sort(fit$subgroups),
    c("CT","MRI")
  )
  
})
