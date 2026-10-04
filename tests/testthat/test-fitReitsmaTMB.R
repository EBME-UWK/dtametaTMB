test_that("fitReitsmaTMB reproduces fitReitsma", {
  
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  
  data("anticcp", package = "dtametaTMB")
  
  fit_glmm <- fitReitsma(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study
  )
  
  fit_tmb <- fitReitsmaTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study
  )
  
  tolerance <- 1e-4
  
  parameter_names <- c(
    "mu_A.sens",
    "mu_B.spec",
    "sigma2_A.sens",
    "sigma2_B.spec",
    "sigma_AB"
  )
  
  # Both optimizers must converge.
  expect_equal(
    fit_glmm$glmmTMB$fit$convergence,
    0
  )
  
  expect_equal(
    fit_tmb$fit$convergence,
    0
  )
  
  # Both Hessians should be positive definite.
  expect_true(
    isTRUE(fit_glmm$glmmTMB$sdr$pdHess)
  )
  
  expect_true(
    isTRUE(fit_tmb$sdreport$pdHess)
  )
  
  # Compare the five fundamental Reitsma parameter estimates.
  expect_equal(
    unname(
      fit_glmm$estimates[
        parameter_names,
        "Estimate"
      ]
    ),
    unname(
      fit_tmb$sdreport2[
        parameter_names,
        "Estimate"
      ]
    ),
    tolerance = tolerance
  )
  
  # Compare delta-method standard errors.
  expect_equal(
    unname(
      fit_glmm$estimates[
        parameter_names,
        "Std_Error"
      ]
    ),
    unname(
      fit_tmb$sdreport2[
        parameter_names,
        "Std. Error"
      ]
    ),
    tolerance = tolerance
  )
  
  # Compare the maximized log-likelihood.
  expect_equal(
    as.numeric(logLik(fit_glmm)),
    as.numeric(logLik(fit_tmb)),
    tolerance = tolerance
  )
  
  # Compare pooled sensitivity and specificity estimates.
  expect_equal(
    unname(fit_glmm$sensspec[, "Estimate"]),
    unname(fit_tmb$sensspec[, "Estimate"]),
    tolerance = tolerance
  )
  
  # Compare sensitivity/specificity confidence limits.
  expect_equal(
    unname(fit_glmm$sensspec[, "CI_Lower"]),
    unname(fit_tmb$sensspec[, "CI_Lower"]),
    tolerance = tolerance
  )
  
  expect_equal(
    unname(fit_glmm$sensspec[, "CI_Upper"]),
    unname(fit_tmb$sensspec[, "CI_Upper"]),
    tolerance = tolerance
  )
  
  # Compare diagnostic odds ratio and likelihood ratios.
  expect_equal(
    unname(fit_glmm$LRDOR[, "Estimate"]),
    unname(fit_tmb$LRDOR[, "Estimate"]),
    tolerance = tolerance
  )
  
  # Compare recovered HSROC parameters.
  hsroc_names <- c(
    "Lambda",
    "Theta",
    "beta",
    "sigma2_alpha",
    "sigma2_theta"
  )
  
  expect_equal(
    unname(
      unlist(
        fit_glmm$RutterGatsonis_recovered[
          1,
          hsroc_names,
          drop = TRUE
        ]
      )
    ),
    unname(
      unlist(
        fit_tmb$RutterGatsonis_recovered[
          1,
          hsroc_names,
          drop = TRUE
        ]
      )
    ),
    tolerance = tolerance
  )
})


test_that("fitReitsmaSubgroupTMB reproduces fitReitsmaSubgroup", {
  
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  
  data("anticcp", package = "dtametaTMB")
  
  fit_glmm <- fitReitsmaSubgroup(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    variances = "common"
  )
  
  fit_tmb <- fitReitsmaSubgroupTMB(
    data = anticcp,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = generation,
    variances = "common"
  )
  
  tolerance <- 1e-4
  
  subgroups <- fit_glmm$subgroups
  subgroup_names_glmm <- make.names(subgroups)
  
  # Both subgroup implementations must preserve the same levels.
  expect_equal(
    fit_tmb$subgroups,
    fit_glmm$subgroups
  )
  
  # Both fitted models must converge.
  expect_equal(
    fit_glmm$glmmTMB_mu$fit$convergence,
    0
  )
  
  expect_equal(
    fit_tmb$fit$convergence,
    0
  )
  
  # Both Hessians should be positive definite.
  expect_true(
    isTRUE(fit_glmm$glmmTMB_mu$sdr$pdHess)
  )
  
  expect_true(
    isTRUE(fit_tmb$sdreport$pdHess)
  )
  
  # Compare subgroup-specific logit sensitivity and specificity.
  for (i in seq_along(subgroups)) {
    
    subgroup_glmm <- subgroup_names_glmm[i]
    subgroup_tmb <- subgroups[i]
    
    mu_A_glmm <- paste0(
      "mu_A.",
      subgroup_glmm
    )
    
    mu_B_glmm <- paste0(
      "mu_B.",
      subgroup_glmm
    )
    
    mu_A_tmb <- paste0(
      "mu_A.",
      subgroup_tmb
    )
    
    mu_B_tmb <- paste0(
      "mu_B.",
      subgroup_tmb
    )
    
    expect_equal(
      unname(
        fit_glmm$estimates_mu[
          c(mu_A_glmm, mu_B_glmm),
          "Estimate"
        ]
      ),
      unname(
        fit_tmb$sdreport2[
          c(mu_A_tmb, mu_B_tmb),
          "Estimate"
        ]
      ),
      tolerance = tolerance,
      label = paste(
        "subgroup means for",
        subgroups[i]
      )
    )
    
    expect_equal(
      unname(
        fit_glmm$estimates_mu[
          c(mu_A_glmm, mu_B_glmm),
          "Std_Error"
        ]
      ),
      unname(
        fit_tmb$sdreport2[
          c(mu_A_tmb, mu_B_tmb),
          "Std. Error"
        ]
      ),
      tolerance = tolerance,
      label = paste(
        "subgroup mean standard errors for",
        subgroups[i]
      )
    )
  }
  
  # Under variances = "common", the glmmTMB implementation reports
  # one common covariance matrix. The custom TMB implementation
  # reports the mapped common values once for each subgroup.
  common_variance_names <- c(
    "sigma2_A.sens",
    "sigma2_B.spec",
    "sigma_AB"
  )
  
  glmm_variance_estimates <- unname(
    fit_glmm$estimates_mu[
      common_variance_names,
      "Estimate"
    ]
  )
  
  glmm_variance_se <- unname(
    fit_glmm$estimates_mu[
      common_variance_names,
      "Std_Error"
    ]
  )
  
  for (subgroup in subgroups) {
    
    tmb_variance_names <- c(
      paste0("sigma2_A.", subgroup),
      paste0("sigma2_B.", subgroup),
      paste0("sigma_AB.", subgroup)
    )
    
    expect_equal(
      glmm_variance_estimates,
      unname(
        fit_tmb$sdreport2[
          tmb_variance_names,
          "Estimate"
        ]
      ),
      tolerance = tolerance,
      label = paste(
        "common variance estimates reported for",
        subgroup
      )
    )
    
    expect_equal(
      glmm_variance_se,
      unname(
        fit_tmb$sdreport2[
          tmb_variance_names,
          "Std. Error"
        ]
      ),
      tolerance = tolerance,
      label = paste(
        "common variance standard errors reported for",
        subgroup
      )
    )
  }
  
  # Compare maximized log-likelihoods.
  expect_equal(
    as.numeric(logLik(fit_glmm)),
    as.numeric(logLik(fit_tmb)),
    tolerance = tolerance
  )
  
  # Align sensitivity/specificity rows before comparing.
  sensspec_glmm <- fit_glmm$sensspec
  sensspec_tmb <- fit_tmb$sensspec
  
  rownames(sensspec_glmm) <- sub(
    paste0(
      "^(",
      paste(subgroup_names_glmm, collapse = "|"),
      ")"
    ),
    "\\1",
    rownames(sensspec_glmm)
  )
  
  expected_glmm_rows <- unlist(
    lapply(
      subgroup_names_glmm,
      function(subgroup) {
        c(
          paste0("mu_A.", subgroup),
          paste0("mu_B.", subgroup)
        )
      }
    ),
    use.names = FALSE
  )
  
  expected_tmb_rows <- unlist(
    lapply(
      subgroups,
      function(subgroup) {
        c(
          paste0("mu_A.", subgroup),
          paste0("mu_B.", subgroup)
        )
      }
    ),
    use.names = FALSE
  )
  
  expect_equal(
    unname(
      sensspec_glmm[
        expected_glmm_rows,
        "Estimate"
      ]
    ),
    unname(
      sensspec_tmb[
        expected_tmb_rows,
        "Estimate"
      ]
    ),
    tolerance = tolerance
  )
  
  expect_equal(
    unname(
      sensspec_glmm[
        expected_glmm_rows,
        "CI_Lower"
      ]
    ),
    unname(
      sensspec_tmb[
        expected_tmb_rows,
        "CI_Lower"
      ]
    ),
    tolerance = tolerance
  )
  
  expect_equal(
    unname(
      sensspec_glmm[
        expected_glmm_rows,
        "CI_Upper"
      ]
    ),
    unname(
      sensspec_tmb[
        expected_tmb_rows,
        "CI_Upper"
      ]
    ),
    tolerance = tolerance
  )
  
  # Compare subgroup-specific recovered HSROC parameters.
  hsroc_names <- c(
    "Lambda",
    "Theta",
    "beta",
    "sigma2_alpha",
    "sigma2_theta"
  )
  
  expect_equal(
    unname(
      as.matrix(
        fit_glmm$RutterGatsonis_recovered[
          subgroups,
          hsroc_names,
          drop = FALSE
        ]
      )
    ),
    unname(
      as.matrix(
        fit_tmb$RutterGatsonis_recovered[
          subgroups,
          hsroc_names,
          drop = FALSE
        ]
      )
    ),
    tolerance = tolerance
  )
})

expect_abs_equal <- function(
    actual,
    expected,
    tolerance,
    label = NULL) {
  
  # Preserve matrix ordering while converting data frames and
  # matrices into ordinary numeric vectors.
  actual <- as.numeric(as.matrix(actual))
  expected <- as.numeric(as.matrix(expected))
  
  testthat::expect_equal(
    length(actual),
    length(expected),
    info = label
  )
  
  testthat::expect_true(
    all(is.finite(actual)),
    info = paste0(
      if (is.null(label)) "" else paste0(label, ": "),
      "actual values contain non-finite values."
    )
  )
  
  testthat::expect_true(
    all(is.finite(expected)),
    info = paste0(
      if (is.null(label)) "" else paste0(label, ": "),
      "expected values contain non-finite values."
    )
  )
  
  differences <- abs(actual - expected)
  max_difference <- max(differences)
  
  testthat::expect_true(
    max_difference <= tolerance,
    info = paste0(
      if (is.null(label)) "" else paste0(label, ": "),
      "maximum absolute difference was ",
      format(
        max_difference,
        scientific = TRUE
      ),
      "; allowed tolerance was ",
      format(
        tolerance,
        scientific = TRUE
      )
    )
  )
  
  invisible(max_difference)
}


extract_glmm_subgroup_parameters <- function(fit) {
  
  subgroups <- fit$subgroups
  subgroup_safe <- make.names(subgroups)
  
  means <- unlist(
    lapply(
      subgroup_safe,
      function(subgroup) {
        fit$estimates_mu[
          c(
            paste0("mu_A.", subgroup),
            paste0("mu_B.", subgroup)
          ),
          "Estimate"
        ]
      }
    ),
    use.names = FALSE
  )
  
  if (fit$variances == "common") {
    
    variances <- fit$estimates_mu[
      c(
        "sigma2_A.sens",
        "sigma2_B.spec",
        "sigma_AB"
      ),
      "Estimate"
    ]
    
  } else {
    
    variances <- unlist(
      lapply(
        subgroup_safe,
        function(subgroup) {
          fit$estimates_mu[
            c(
              paste0("sigma2_A.", subgroup),
              paste0("sigma2_B.", subgroup),
              paste0("sigma_AB.", subgroup)
            ),
            "Estimate"
          ]
        }
      ),
      use.names = FALSE
    )
  }
  
  unname(
    c(
      means,
      variances
    )
  )
}


extract_tmb_subgroup_parameters <- function(fit) {
  
  subgroups <- fit$subgroups
  
  means <- unlist(
    lapply(
      subgroups,
      function(subgroup) {
        fit$sdreport2[
          c(
            paste0("mu_A.", subgroup),
            paste0("mu_B.", subgroup)
          ),
          "Estimate"
        ]
      }
    ),
    use.names = FALSE
  )
  
  if (fit$variances == "common") {
    
    first_subgroup <- subgroups[1]
    
    variances <- fit$sdreport2[
      c(
        paste0("sigma2_A.", first_subgroup),
        paste0("sigma2_B.", first_subgroup),
        paste0("sigma_AB.", first_subgroup)
      ),
      "Estimate"
    ]
    
  } else {
    
    variances <- unlist(
      lapply(
        subgroups,
        function(subgroup) {
          fit$sdreport2[
            c(
              paste0("sigma2_A.", subgroup),
              paste0("sigma2_B.", subgroup),
              paste0("sigma_AB.", subgroup)
            ),
            "Estimate"
          ]
        }
      ),
      use.names = FALSE
    )
  }
  
  unname(
    c(
      means,
      variances
    )
  )
}

test_that(
  "fitReitsmaSubgroupTMB reproduces unequal-variance glmmTMB model",
  {
    
    skip_on_cran()
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("TMB")
    
    data(
      "anticcp",
      package = "dtametaTMB"
    )
    
    fit_glmm <- fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "unequal"
    )
    
    fit_tmb <- fitReitsmaSubgroupTMB(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "unequal"
    )
    
    expect_equal(
      fit_glmm$subgroups,
      fit_tmb$subgroups
    )
    
    expect_equal(
      fit_glmm$glmmTMB_mu$fit$convergence,
      0
    )
    
    expect_equal(
      fit_tmb$fit$convergence,
      0
    )
    
    expect_true(
      isTRUE(fit_glmm$glmmTMB_mu$sdr$pdHess)
    )
    
    expect_true(
      isTRUE(fit_tmb$sdreport$pdHess)
    )
    
    glmm_parameters <-
      extract_glmm_subgroup_parameters(fit_glmm)
    
    tmb_parameters <-
      extract_tmb_subgroup_parameters(fit_tmb)
    
    expect_abs_equal(
      actual = glmm_parameters,
      expected = tmb_parameters,
      tolerance = 1e-4,
      label = paste(
        "unequal-variance Reitsma subgroup parameters",
        "for Anti-CCP"
      )
    )
    
    expect_abs_equal(
      actual = logLik(fit_glmm),
      expected = logLik(fit_tmb),
      tolerance = 1e-5,
      label = paste(
        "unequal-variance log-likelihood",
        "for Anti-CCP"
      )
    )
    
    hsroc_names <- c(
      "Lambda",
      "Theta",
      "beta",
      "sigma2_alpha",
      "sigma2_theta"
    )
    
    expect_abs_equal(
      actual =
        fit_glmm$RutterGatsonis_recovered[
          fit_glmm$subgroups,
          hsroc_names,
          drop = FALSE
        ],
      expected =
        fit_tmb$RutterGatsonis_recovered[
          fit_tmb$subgroups,
          hsroc_names,
          drop = FALSE
        ],
      tolerance = 2e-4,
      label = paste(
        "unequal-variance recovered HSROC parameters",
        "for Anti-CCP"
      )
    )
  }
)


test_that(
  "fitReitsmaSubgroupTMB reproduces common sensitivity model",
  {
    
    skip_on_cran()
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("TMB")
    
    data(
      "anticcp",
      package = "dtametaTMB"
    )
    
    fit_glmm <- fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = "sens"
    )
    
    fit_tmb <- fitReitsmaSubgroupTMB(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = "sens"
    )
    
    expect_equal(
      fit_glmm$glmmTMB_mu$fit$convergence,
      0
    )
    
    expect_equal(
      fit_tmb$fit$convergence,
      0
    )
    
    expect_true(
      isTRUE(fit_glmm$glmmTMB_mu$sdr$pdHess)
    )
    
    expect_true(
      isTRUE(fit_tmb$sdreport$pdHess)
    )
    
    glmm_parameters <-
      extract_glmm_subgroup_parameters(fit_glmm)
    
    tmb_parameters <-
      extract_tmb_subgroup_parameters(fit_tmb)
    
    expect_abs_equal(
      actual = glmm_parameters,
      expected = tmb_parameters,
      tolerance = 1e-4,
      label = "common-sensitivity Reitsma subgroup parameters"
    )
    
    expect_abs_equal(
      actual = logLik(fit_glmm),
      expected = logLik(fit_tmb),
      tolerance = 1e-5,
      label = "common-sensitivity log-likelihood"
    )
    
    glmm_sensitivity <- fit_glmm$estimates_mu[
      grep(
        "^mu_A\\.",
        rownames(fit_glmm$estimates_mu)
      ),
      "Estimate"
    ]
    
    tmb_sensitivity <- fit_tmb$sdreport2[
      grep(
        "^mu_A\\.",
        rownames(fit_tmb$sdreport2)
      ),
      "Estimate"
    ]
    
    expect_lt(
      diff(range(glmm_sensitivity)),
      1e-10
    )
    
    expect_lt(
      diff(range(tmb_sensitivity)),
      1e-10
    )
  }
)


test_that(
  "fitReitsmaSubgroupTMB reproduces common specificity model",
  {
    
    skip_on_cran()
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("TMB")
    
    data(
      "anticcp",
      package = "dtametaTMB"
    )
    
    fit_glmm <- fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = "spec"
    )
    
    fit_tmb <- fitReitsmaSubgroupTMB(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = "spec"
    )
    
    expect_equal(
      fit_glmm$glmmTMB_mu$fit$convergence,
      0
    )
    
    expect_equal(
      fit_tmb$fit$convergence,
      0
    )
    
    expect_true(
      isTRUE(fit_glmm$glmmTMB_mu$sdr$pdHess)
    )
    
    expect_true(
      isTRUE(fit_tmb$sdreport$pdHess)
    )
    
    expect_abs_equal(
      actual =
        extract_glmm_subgroup_parameters(fit_glmm),
      expected =
        extract_tmb_subgroup_parameters(fit_tmb),
      tolerance = 1e-4,
      label = "common-specificity Reitsma subgroup parameters"
    )
    
    expect_abs_equal(
      actual = logLik(fit_glmm),
      expected = logLik(fit_tmb),
      tolerance = 1e-5,
      label = "common-specificity log-likelihood"
    )
    
    glmm_specificity <- fit_glmm$estimates_mu[
      grep(
        "^mu_B\\.",
        rownames(fit_glmm$estimates_mu)
      ),
      "Estimate"
    ]
    
    tmb_specificity <- fit_tmb$sdreport2[
      grep(
        "^mu_B\\.",
        rownames(fit_tmb$sdreport2)
      ),
      "Estimate"
    ]
    
    expect_lt(
      diff(range(glmm_specificity)),
      1e-10
    )
    
    expect_lt(
      diff(range(tmb_specificity)),
      1e-10
    )
  }
)

test_that(
  "fitReitsmaSubgroupTMB reproduces model with common sensitivity and specificity",
  {
    
    skip_on_cran()
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("TMB")
    
    data(
      "anticcp",
      package = "dtametaTMB"
    )
    
    constraints <- c(
      "sens",
      "spec"
    )
    
    fit_glmm <- fitReitsmaSubgroup(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = constraints
    )
    
    fit_tmb <- fitReitsmaSubgroupTMB(
      data = anticcp,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      study = study,
      subgroup = generation,
      variances = "common",
      sensspec_constrain = constraints
    )
    
    expect_equal(
      fit_glmm$glmmTMB_mu$fit$convergence,
      0
    )
    
    expect_equal(
      fit_tmb$fit$convergence,
      0
    )
    
    expect_true(
      isTRUE(fit_glmm$glmmTMB_mu$sdr$pdHess)
    )
    
    expect_true(
      isTRUE(fit_tmb$sdreport$pdHess)
    )
    
    expect_abs_equal(
      actual =
        extract_glmm_subgroup_parameters(fit_glmm),
      expected =
        extract_tmb_subgroup_parameters(fit_tmb),
      tolerance = 1e-4,
      label = paste(
        "common-sensitivity-and-specificity",
        "Reitsma subgroup parameters"
      )
    )
    
    expect_abs_equal(
      actual = logLik(fit_glmm),
      expected = logLik(fit_tmb),
      tolerance = 1e-5,
      label = paste(
        "common-sensitivity-and-specificity",
        "log-likelihood"
      )
    )
    
    glmm_sensitivity <- fit_glmm$estimates_mu[
      grep(
        "^mu_A\\.",
        rownames(fit_glmm$estimates_mu)
      ),
      "Estimate"
    ]
    
    glmm_specificity <- fit_glmm$estimates_mu[
      grep(
        "^mu_B\\.",
        rownames(fit_glmm$estimates_mu)
      ),
      "Estimate"
    ]
    
    tmb_sensitivity <- fit_tmb$sdreport2[
      grep(
        "^mu_A\\.",
        rownames(fit_tmb$sdreport2)
      ),
      "Estimate"
    ]
    
    tmb_specificity <- fit_tmb$sdreport2[
      grep(
        "^mu_B\\.",
        rownames(fit_tmb$sdreport2)
      ),
      "Estimate"
    ]
    
    expect_lt(
      diff(range(glmm_sensitivity)),
      1e-10
    )
    
    expect_lt(
      diff(range(glmm_specificity)),
      1e-10
    )
    
    expect_lt(
      diff(range(tmb_sensitivity)),
      1e-10
    )
    
    expect_lt(
      diff(range(tmb_specificity)),
      1e-10
    )
  }
)


test_that(
  "both Reitsma subgroup backends reject variance constraints with unequal variances",
  {
    
    data(
      "anticcp",
      package = "dtametaTMB"
    )
    
    expect_error(
      fitReitsmaSubgroup(
        data = anticcp,
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        study = study,
        subgroup = generation,
        variances = "unequal",
        constrain = "sigma_AB"
      ),
      "currently not supported"
    )
    
    expect_error(
      fitReitsmaSubgroupTMB(
        data = anticcp,
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        study = study,
        subgroup = generation,
        variances = "unequal",
        constrain = "sigma_AB"
      ),
      "currently not supported"
    )
  }
)


