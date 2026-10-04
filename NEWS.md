# dtametaTMB 0.1.4

* Added `symbol` argument to all `plot.*()` methods, allowing study-level
  markers to be drawn as rectangles, ellipses, diamonds, triangles,
  crosses, plus-signs, or stars.
* Added `size = "sampsize_revman"` and `size = "se_revman"`,
  `size = eb_revman` options, following RevMan-style symbol-sizing convention, 
  to support visual consistency with `as_revman()` exports.
* Corrected bug `size = "fisher"` which the correct covariance now to construct 
  the matrix G for the SROC plots.
* Corrected bug in `plot` and `as_revman` methods for ReitsmaSubgroup models when 
  sensitivity and/or specificity are assumed common in subgroups.
* Fixed manual page titles to consistently use title case.

# dtametaTMB 0.1.3

* `summary()` methods now return objects of a dedicated S3 class
  (`"summary.DTAmodel"` and model-specific subclasses) rather than plain
  lists, with a corresponding `print.summary.DTAmodel()` method.
* Refactored forest plot construction across model classes into shared
  internal helpers (`getForestSensSpec()`, `getForestSensSpecLCA()`,
  `getForestPlot()`) to reduce code duplication.
* `coef()` and `vcov()` methods are now available for all model classes;
  examples updated to use `summary()` and these accessors instead of
  accessing internal list components directly.

# dtametaTMB 0.1.2

* Maintenance release.

# dtametaTMB 0.1.1

* Added latent class models.

# dtametaTMB 0.1.0

* Initial CRAN release.