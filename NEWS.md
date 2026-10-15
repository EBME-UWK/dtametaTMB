# dtametaTMB 0.1.4

* Added a `symbol` argument to all `plot.*()` methods, allowing study-level
  markers to be drawn as rectangles, ellipses, diamonds, triangles, crosses,
  plus signs, or stars.

* Subgroup plots now distinguish subgroups by point symbol as well as colour
  by default. The previous appearance, using rectangles for all subgroups,
  can be reproduced by setting `symbol = "rectangle"`.

* Added the `size = "sampsize_revman"` and `size = "se_revman"` options,
  following RevMan-style symbol-sizing conventions, to support visual
  consistency with `as_revman()` exports. The analogous
  `size = "eb_revman"` option is available for latent-class models.

* Corrected a bug in `plot.Reitsma()` with the default `size = "fisher"`.
  The covariance used to construct the matrix `G` for the SROC plot is now
  calculated correctly.

* Added the `size = "fisher"` option to `plot.RutterGatsonis()`,
  `plot.RutterGatsonisSubgroup()`, and `plot.ReitsmaSubgroup()`. It is now
  the default in `plot.RutterGatsonis()`.

* Added `fitReitsmaTMB()` and `fitReitsmaSubgroupTMB()`, which use custom
  TMB templates rather than `glmmTMB` and may provide improved numerical
  robustness in difficult or near-boundary cases.

* Corrected a bug in the `plot()` and `as_revman()` methods for
  `ReitsmaSubgroup` models when sensitivity and/or specificity are
  constrained to be common across subgroups.

* Added optional study-level confidence intervals (`studyCI = TRUE`) and
  study labels (`studylabels = TRUE`) to SROC plotting methods.

* Improved starting-value initialization across model families. Robust,
  finite fallback values are now used when preliminary model fits or
  empirical estimates are unavailable, including for sparse and degenerate
  datasets.

* Revised the likelihood-ratio-test examples and documentation to caution
  against using the ordinary chi-squared reference distribution when a
  reduced model constrains one or more random-effects variances to zero.
  The corresponding zero-variance boundary comparison was removed from the
  vignette.

* Fixed manual-page titles to use title case consistently.

# dtametaTMB 0.1.3

* `summary()` methods now return objects of dedicated S3 classes
  (`"summary.DTAmodel"` and model-specific subclasses) rather than plain
  lists, with a corresponding `print.summary.DTAmodel()` method.

* Refactored forest-plot construction across model classes into shared
  internal helpers (`getForestSensSpec()`, `getForestSensSpecLCA()`, and
  `getForestPlot()`) to reduce code duplication.

* `coef()` and `vcov()` methods are now available for all model classes.
  Examples were updated to use `summary()` and these accessors rather than
  accessing internal list components directly.

# dtametaTMB 0.1.2

* Maintenance release.

# dtametaTMB 0.1.1

* Added latent-class models.

# dtametaTMB 0.1.0

* Initial CRAN release.