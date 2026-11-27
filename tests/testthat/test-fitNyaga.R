# -------------------------------------------------------------------------
# 0. Load package datasets and prepare standardised test data
# -------------------------------------------------------------------------
skip_on_cran()

.nyaga_dataset_env <- new.env(parent = emptyenv())
utils::data("HPV", package = "dtametaTMB", envir = .nyaga_dataset_env)
utils::data("anticcp", package = "dtametaTMB", envir = .nyaga_dataset_env)

# Keep the reference test first, matching the NMA vignette output.
.nyaga_hpv_raw <- .nyaga_dataset_env$HPV
.nyaga_hpv <- data.frame(
  study = .nyaga_hpv_raw$StudyID,
  test = factor(.nyaga_hpv_raw$Test, levels = c("HC2", "CC", "LBC")),
  TP = .nyaga_hpv_raw$TP,
  FP = .nyaga_hpv_raw$FP,
  FN = .nyaga_hpv_raw$FN,
  TN = .nyaga_hpv_raw$TN
)

# anticcp contains one study per row. Generation (CCP1/CCP2) is used as
# the two-level test/group factor for the model-reduction identity.
.nyaga_anticcp_raw <- .nyaga_dataset_env$anticcp
.nyaga_anticcp <- data.frame(
  study = .nyaga_anticcp_raw$study,
  test = factor(.nyaga_anticcp_raw$generation,
                levels = c("CCP1", "CCP2")),
  TP = .nyaga_anticcp_raw$TP,
  FP = .nyaga_anticcp_raw$FP,
  FN = .nyaga_anticcp_raw$FN,
  TN = .nyaga_anticcp_raw$TN
)

# fitNyaga() uses substitute() to capture column names, so these wrappers
# deliberately pass the column names as unquoted symbols.
.fit_nyaga_test <- function(dat, ...) {
  fitNyaga(
    data = dat,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    test = test,
    ...
  )
}

.fit_reitsma_subgroup_test <- function(dat, ...) {
  fitReitsmaSubgroupTMB(
    data = dat,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    subgroup = test,
    ...
  )
}

# Cache fits used in multiple tests. Each key corresponds to exactly one
# model specification; this avoids repeatedly fitting identical TMB models.
.nyaga_fit_cache <- new.env(parent = emptyenv())

.get_nyaga_fit <- function(key, dat = .nyaga_hpv, ...) {
  if (!exists(key, envir = .nyaga_fit_cache, inherits = FALSE)) {
    assign(key, .fit_nyaga_test(dat, ...), envir = .nyaga_fit_cache)
  }
  get(key, envir = .nyaga_fit_cache, inherits = FALSE)
}

# Run graphics tests on a temporary PDF device, including on headless systems.
.run_nyaga_pdf <- function(expr) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  force(expr)
}

# -------------------------------------------------------------------------
# 1. Check the bundled datasets used by the tests
# -------------------------------------------------------------------------

test_that("bundled HPV and anticcp datasets have the expected structures", {
  skip_on_cran()
  expect_equal(nrow(.nyaga_hpv), 36L)
  expect_equal(length(unique(.nyaga_hpv$study)), 17L)
  expect_equal(levels(.nyaga_hpv$test), c("HC2", "CC", "LBC"))
  expect_true(all(table(.nyaga_hpv$study) >= 2L))
  expect_equal(
    sum(duplicated(.nyaga_hpv[c("study", "test")])),
    0L
  )

  expect_equal(nrow(.nyaga_anticcp), 37L)
  expect_equal(length(unique(.nyaga_anticcp$study)), 37L)
  expect_equal(levels(.nyaga_anticcp$test), c("CCP1", "CCP2"))
  expect_true(all(table(.nyaga_anticcp$study) == 1L))
})

# -------------------------------------------------------------------------
# 2. Input validation
# -------------------------------------------------------------------------

test_that("fitNyaga requires at least two diagnostic tests", {
  skip_on_cran()
  one_test <- .nyaga_hpv[.nyaga_hpv$test == "HC2", , drop = FALSE]

  expect_error(
    .fit_nyaga_test(one_test),
    "requires at least two diagnostic tests"
  )
})

test_that("fitNyaga rejects duplicate study-test combinations", {
  skip_on_cran()
  duplicated_arm <- rbind(
    .nyaga_hpv,
    .nyaga_hpv[1L, , drop = FALSE]
  )

  expect_error(
    .fit_nyaga_test(duplicated_arm),
    "Each study-test combination must occur only once"
  )
})

test_that("fitNyaga validates diagnostic accuracy counts", {
  skip_on_cran()
  bad_data <- .nyaga_hpv
  bad_data$TP[1L] <- -1L
  expect_error(.fit_nyaga_test(bad_data))

})

test_that("fitNyaga validates confidence levels", {
  skip_on_cran()
  expect_error(.fit_nyaga_test(.nyaga_hpv, conflevel = 0))
  expect_error(.fit_nyaga_test(.nyaga_hpv, conflevel = 1))
  expect_error(.fit_nyaga_test(.nyaga_hpv, conflevel = c(0.90, 0.95)))
})

test_that("fitNyaga validates study-level, arm-level and fixed-effect constraints", {
  skip_on_cran()
  expect_error(.fit_nyaga_test(.nyaga_hpv, constrain = "unknown"))
  expect_error(.fit_nyaga_test(.nyaga_hpv, constrain = c("sigma_AB", "all")))
  expect_error(.fit_nyaga_test(.nyaga_hpv, tau_constrain = "tau_A"))
  expect_error(.fit_nyaga_test(.nyaga_hpv, tau_constrain = 1))
  expect_error(.fit_nyaga_test(.nyaga_hpv, sensspec_constrain = "accuracy"))
  expect_error(.fit_nyaga_test(.nyaga_hpv, sensspec_constrain = 1))
})

# -------------------------------------------------------------------------
# 3. Main multi-arm HPV model: convergence and regression benchmarks
# -------------------------------------------------------------------------

test_that("fitNyaga fits the bundled HPV network and returns the expected structure", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")

  expect_s3_class(fit, "Nyaga")
  expect_s3_class(fit, "DTAmodel")
  expect_equal(fit$tests, c("HC2", "CC", "LBC"))
  expect_equal(nrow(fit$data), 36L)
  expect_equal(nlevels(fit$data$study), 17L)

  expect_equal(as.integer(fit$fit$convergence), 0L)
  expect_true(isTRUE(fit$sdreport$pdHess))
  expect_lt(max(abs(fit$sdreport$gradient.fixed)), 1e-3)

  required_components <- c(
    "data", "fit", "sdreport", "sdreport2", "vcov",
    "sensspec", "LRDOR", "tests", "constrain",
    "tau_constrain", "sensspec_constrain"
  )
  expect_true(all(required_components %in% names(fit)))

  expected_parameter_rows <- c(
    paste0("mu_A.", fit$tests),
    paste0("mu_B.", fit$tests),
    "sigma2_A", "sigma2_B", "sigma_AB", "rho_AB",
    "tau2_A", "tau2_B",
    "total_var_A", "total_var_B",
    "icc_A", "icc_B", "total_rho_AB",
    paste0("nu_A.", fit$tests[-1L]),
    paste0("nu_B.", fit$tests[-1L])
  )

  expect_true(all(expected_parameter_rows %in% rownames(fit$sdreport2)))
  expect_true(all(is.finite(fit$sdreport2[, "Estimate"])))
  expect_true(all(is.finite(fit$sdreport2[, "Std. Error"])))
})

test_that("fitNyaga reproduces the HPV vignette parameter benchmarks", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")
  cf <- stats::coef(fit)

  # Reference estimates are from the NMA vignette built with this package
  # version. Tolerance allows for small numerical differences across systems.
  expected <- c(
    "mu_A.HC2" = 2.7564231,
    "mu_A.CC" = 0.54896379,
    "mu_A.LBC" = 1.4350793,
    "mu_B.HC2" = 2.1224412,
    "mu_B.CC" = 3.0136616,
    "mu_B.LBC" = 2.4537805,
    "sigma2_A" = 0.082036906,
    "sigma2_B" = 0.19397243,
    "sigma_AB" = -0.070914053,
    "tau2_A" = 0.66967591,
    "tau2_B" = 0.19852073
  )

  expect_true(all(names(expected) %in% names(cf)))
  expect_equal(
    unname(cf[names(expected)]),
    unname(expected),
    tolerance = 5e-3
  )
  expect_equal(
    as.numeric(stats::logLik(fit)),
    -357.5425,
    tolerance = 5e-3
  )
  expect_equal(attr(stats::logLik(fit), "df"), 11L)
  expect_equal(attr(stats::logLik(fit), "nobs"), 36L)
})

test_that("fitNyaga returns coherent sensitivity, specificity and LRDOR summaries", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")

  expect_true(is.data.frame(fit$sensspec))
  expect_true(is.data.frame(fit$LRDOR))
  expect_equal(nrow(fit$sensspec), 2L * length(fit$tests))
  expect_equal(nrow(fit$LRDOR), 3L * length(fit$tests))

  expect_true(all(c(
    "type", "Estimate", "conflevel", "CI_Lower", "CI_Upper"
  ) %in% names(fit$sensspec)))

  expect_true(all(fit$sensspec$Estimate > 0 & fit$sensspec$Estimate < 1))
  expect_true(all(fit$sensspec$CI_Lower <= fit$sensspec$Estimate))
  expect_true(all(fit$sensspec$Estimate <= fit$sensspec$CI_Upper))
  expect_true(all(fit$sensspec$CI_Lower >= 0 & fit$sensspec$CI_Upper <= 1))

  expected_lrdor_rows <- unlist(lapply(
    fit$tests,
    function(z) paste0(z, ": ", c("DOR", "LR+", "LR-"))
  ))
  expect_equal(rownames(fit$LRDOR), expected_lrdor_rows)
  expect_true(all(fit$LRDOR$Estimate > 0))
  expect_true(all(fit$LRDOR$CI_Lower > 0))
  expect_true(all(fit$LRDOR$CI_Lower <= fit$LRDOR$Estimate))
  expect_true(all(fit$LRDOR$Estimate <= fit$LRDOR$CI_Upper))

  for (test_name in fit$tests) {
    mu_A <- fit$sdreport2[paste0("mu_A.", test_name), "Estimate"]
    mu_B <- fit$sdreport2[paste0("mu_B.", test_name), "Estimate"]
    sens <- stats::plogis(mu_A)
    spec <- stats::plogis(mu_B)

    expect_equal(
      fit$LRDOR[paste0(test_name, ": DOR"), "Estimate"],
      exp(mu_A + mu_B),
      tolerance = 1e-7
    )
    expect_equal(
      fit$LRDOR[paste0(test_name, ": LR+"), "Estimate"],
      sens / (1 - spec),
      tolerance = 1e-7
    )
    expect_equal(
      fit$LRDOR[paste0(test_name, ": LR-"), "Estimate"],
      (1 - sens) / spec,
      tolerance = 1e-7
    )
  }
})

test_that("fitNyaga reports internally consistent heterogeneity quantities", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")
  b <- fit$sdreport2
  est <- function(parameter) unname(b[parameter, "Estimate"])

  expect_equal(est("total_var_A"), est("sigma2_A") + est("tau2_A"),
               tolerance = 1e-7)
  expect_equal(est("total_var_B"), est("sigma2_B") + est("tau2_B"),
               tolerance = 1e-7)
  expect_equal(est("icc_A"), est("sigma2_A") / est("total_var_A"),
               tolerance = 1e-7)
  expect_equal(est("icc_B"), est("sigma2_B") / est("total_var_B"),
               tolerance = 1e-7)
  expect_equal(
    est("total_rho_AB"),
    est("sigma_AB") / sqrt(est("total_var_A") * est("total_var_B")),
    tolerance = 1e-7
  )

  for (test_name in fit$tests[-1L]) {
    expect_equal(
      est(paste0("nu_A.", test_name)),
      est(paste0("mu_A.", test_name)) - est(paste0("mu_A.", fit$tests[1L])),
      tolerance = 1e-7
    )
    expect_equal(
      est(paste0("nu_B.", test_name)),
      est(paste0("mu_B.", test_name)) - est(paste0("mu_B.", fit$tests[1L])),
      tolerance = 1e-7
    )
  }
})

# -------------------------------------------------------------------------
# 4. S3 methods and reference-test invariance
# -------------------------------------------------------------------------

test_that("Nyaga coefficients, covariance matrix and likelihood methods are consistent", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")

  cf <- stats::coef(fit)
  vc <- stats::vcov(fit)
  ll <- stats::logLik(fit)
  sm <- summary(fit)

  expect_equal(unname(cf), unname(fit$sdreport2[, "Estimate"]),
               tolerance = 1e-10)
  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), length(cf))
  expect_equal(ncol(vc), length(cf))
  expect_equal(rownames(vc), names(cf))
  expect_equal(colnames(vc), names(cf))
  expect_equal(vc, t(vc), tolerance = 1e-7)

  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -fit$fit$objective, tolerance = 1e-10)
  expect_equal(attr(ll, "df"), length(fit$fit$par))
  expect_equal(attr(ll, "nobs"), nrow(fit$data))
  expect_true(is.finite(stats::AIC(fit)))
  expect_true(is.finite(stats::BIC(fit)))

  expect_true(all(c("summary.Nyaga", "summary.DTAmodel") %in% class(sm)))
  expect_named(sm, c("estimates", "sensspec", "LRDOR","RutterGatsonis_recovered", "tests"))
  expect_equal(sm$tests, fit$tests)
  expect_equal(sm$sensspec, fit$sensspec)
  expect_equal(sm$LRDOR, fit$LRDOR)
  expect_output(print(fit), "Nyaga NMA Model")
  expect_no_error(print(sm))
})

test_that("changing the reference test leaves fitted test-specific means unchanged", {
  skip_on_cran()
  fit_hc2_reference <- .get_nyaga_fit("hpv_default")

  dat_cc_reference <- .nyaga_hpv
  dat_cc_reference$test <- factor(
    dat_cc_reference$test,
    levels = c("CC", "HC2", "LBC")
  )
  fit_cc_reference <- .get_nyaga_fit(
    "hpv_cc_reference",
    dat = dat_cc_reference
  )

  cf_hc2 <- stats::coef(fit_hc2_reference)
  cf_cc <- stats::coef(fit_cc_reference)
  mean_names <- c(
    "mu_A.HC2", "mu_A.CC", "mu_A.LBC",
    "mu_B.HC2", "mu_B.CC", "mu_B.LBC"
  )

  expect_equal(unname(cf_hc2[mean_names]), unname(cf_cc[mean_names]),
               tolerance = 5e-4)
  expect_equal(as.numeric(stats::logLik(fit_hc2_reference)),
               as.numeric(stats::logLik(fit_cc_reference)),
               tolerance = 5e-4)
  expect_equal(fit_cc_reference$tests, c("CC", "HC2", "LBC"))
})

# -------------------------------------------------------------------------
# 5. Constraint handling
# -------------------------------------------------------------------------

test_that("Nyaga study-level covariance and variance constraints are applied", {
  skip_on_cran()
  fit_cov0 <- .get_nyaga_fit("hpv_sigma_AB_zero", constrain = "sigma_AB")
  expect_equal(unname(stats::coef(fit_cov0)["sigma_AB"]), 0,
               tolerance = 1e-12)

  fit_varA0 <- .get_nyaga_fit("hpv_sigma2_A_zero", constrain = "sigma2_A")
  expect_lt(abs(unname(stats::coef(fit_varA0)["sigma2_A"])), 1e-24)
  expect_equal(unname(stats::coef(fit_varA0)["sigma_AB"]), 0,
               tolerance = 1e-12)

  fit_varB0 <- .get_nyaga_fit("hpv_sigma2_B_zero", constrain = "sigma2_B")
  expect_lt(abs(unname(stats::coef(fit_varB0)["sigma2_B"])), 1e-24)
  expect_equal(unname(stats::coef(fit_varB0)["sigma_AB"]), 0,
               tolerance = 1e-12)

  fit_all0 <- .get_nyaga_fit("hpv_all_study_variances_zero", constrain = "all")
  expect_lt(abs(unname(stats::coef(fit_all0)["sigma2_A"])), 1e-24)
  expect_lt(abs(unname(stats::coef(fit_all0)["sigma2_B"])), 1e-24)
  expect_equal(unname(stats::coef(fit_all0)["sigma_AB"]), 0,
               tolerance = 1e-12)
})

test_that("Nyaga arm-level variance constraints are applied", {
  # Use the real multi-arm HPV network. The two constraints are fitted
  # separately so each mapping is exercised independently.
  skip_on_cran()
  fit_tauA0 <- .get_nyaga_fit("hpv_tau2_A_zero", tau_constrain = "tau2_A")
  expect_lt(abs(unname(stats::coef(fit_tauA0)["tau2_A"])), 1e-24)

  fit_tauB0 <- .get_nyaga_fit("hpv_tau2_B_zero", tau_constrain = "tau2_B")
  expect_lt(abs(unname(stats::coef(fit_tauB0)["tau2_B"])), 1e-24)
})

test_that("Nyaga sensitivity and specificity constraints are applied", {
  skip_on_cran()
  fit_equal_sens <- .get_nyaga_fit(
    "hpv_equal_sensitivity",
    sensspec_constrain = "sens"
  )
  fit_equal_spec <- .get_nyaga_fit(
    "hpv_equal_specificity",
    sensspec_constrain = "spec"
  )
  fit_equal_both <- .get_nyaga_fit(
    "hpv_equal_sensitivity_specificity",
    sensspec_constrain = c("sens", "spec")
  )

  check_equal_parameters <- function(fit, parameter_names) {
    cf <- stats::coef(fit)
    expect_equal(
      unname(cf[parameter_names]),
      rep(unname(cf[parameter_names[1L]]), length(parameter_names)),
      tolerance = 1e-8
    )
  }

  check_equal_parameters(
    fit_equal_sens,
    paste0("mu_A.", fit_equal_sens$tests)
  )
  check_equal_parameters(
    fit_equal_spec,
    paste0("mu_B.", fit_equal_spec$tests)
  )
  check_equal_parameters(
    fit_equal_both,
    paste0("mu_A.", fit_equal_both$tests)
  )
  check_equal_parameters(
    fit_equal_both,
    paste0("mu_B.", fit_equal_both$tests)
  )
})

test_that("Nyaga fixed-effect restrictions can be compared by likelihood-ratio test", {
  skip_on_cran()
  fit_full <- .get_nyaga_fit("hpv_default")
  fit_restricted <- .get_nyaga_fit(
    "hpv_equal_sensitivity_specificity",
    sensspec_constrain = c("sens", "spec")
  )

  comparison <- stats::anova(fit_restricted, fit_full)

  expect_s3_class(comparison, "anova")
  expect_equal(nrow(comparison), 2L)
  expect_gt(comparison$Df.diff[2L], 0)
  expect_true(is.finite(comparison$Chisq[2L]))
})

# -------------------------------------------------------------------------
# 6. Model-reduction identity using anticcp
# -------------------------------------------------------------------------

test_that(paste(
  "Nyaga reduces to the common-variance Reitsma subgroup model",
  "for one test per study when both tau variances are fixed to zero"
), {
  skip_on_cran()
  fit_nyaga_reduced <- .fit_nyaga_test(
    .nyaga_anticcp,
    tau_constrain = c("tau2_A", "tau2_B")
  )

  fit_reitsma_subgroup <- .fit_reitsma_subgroup_test(
    .nyaga_anticcp,
    variances = "common"
  )

  test_levels <- levels(.nyaga_anticcp$test)
  mean_names <- c(
    paste0("mu_A.", test_levels),
    paste0("mu_B.", test_levels)
  )

  expect_equal(
    unname(stats::coef(fit_nyaga_reduced)[mean_names]),
    unname(stats::coef(fit_reitsma_subgroup)[mean_names]),
    tolerance = 2e-3
  )

  # In the common-variance subgroup fit, the shared covariance estimates
  # are repeated under the subgroup-specific report names.
  for (parameter in c("sigma2_A", "sigma2_B", "sigma_AB")) {
    subgroup_names <- paste0(parameter, ".", test_levels)
    nyaga_value <- unname(stats::coef(fit_nyaga_reduced)[parameter])
    subgroup_values <- unname(stats::coef(fit_reitsma_subgroup)[subgroup_names])

    expect_equal(
      rep(nyaga_value, length(test_levels)),
      subgroup_values,
      tolerance = 2e-3
    )
  }

  expect_lt(abs(unname(stats::coef(fit_nyaga_reduced)["tau2_A"])), 1e-24)
  expect_lt(abs(unname(stats::coef(fit_nyaga_reduced)["tau2_B"])), 1e-24)

  expect_equal(
    as.numeric(stats::logLik(fit_nyaga_reduced)),
    as.numeric(stats::logLik(fit_reitsma_subgroup)),
    tolerance = 5e-3
  )
})

# -------------------------------------------------------------------------
# 7. Graphical methods
# -------------------------------------------------------------------------

test_that("Nyaga summary ROC and forest plots run on the HPV network", {
  skip_on_cran()
  fit <- .get_nyaga_fit("hpv_default")

  expect_no_error(.run_nyaga_pdf(
    plot(
      fit,
      size = "sampsize_revman",
      studyCI = TRUE,
      studylabels = TRUE,
      legend = FALSE
    )
  ))

  forest_object <- NULL
  expect_no_error(
    forest_object <- .run_nyaga_pdf(forest(fit, order = "test"))
  )
  expect_false(is.null(forest_object))

  expect_error(plot(fit, conflevel = 1), "conflevel must be a single number")
  expect_error(plot(fit, predlevel = 0), "predlevel must be a single number")
  expect_error(forest(fit, conflevel = 1), "conflevel must be a single number")
})
