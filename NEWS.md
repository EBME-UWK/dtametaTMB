# dtametaTMB 0.1.4

* Added `symbol` argument to all `plot.*()` methods, allowing study-level
  markers to be drawn as rectangles, ellipses, diamonds, triangles,
  crosses, plus-signs, or stars.
* Subgroup plots now distinguish subgroups by point symbol as well as colour
  by default. The previous appearance, using rectangles for all subgroups, 
  can be reproduced by setting `symbol = "rectangle"`.
* Added `size = "sampsize_revman"` and `size = "se_revman"`,
  `size = "eb_revman"` options, following RevMan-style symbol-sizing convention, 
  to support visual consistency with `as_revman()` exports.
* Corrected a bug in `plot.Reitsma()` with default `size = "fisher"`. 
  The covariance used to construct the matrix G for the SROC plot is now calculated correctly.
* Added `size = "fisher"` option to `plot.RutterGatsonis()`, `plot.RutterGatsonisSubgroup()`,
  and `plot.ReitsmaSubgroup()`, which is new default in `plot.RutterGatsonis()` and `plot.ReitsmaSubgroup()`.
* Added `fitReitsmaTMB()` and `fitReitsmaSubgroupTMB()`, which rely on custom TMB templates rather
  than on `glmmTMB` and may offer improved numerical robustness in difficult or near-boundary cases.
* Corrected a bug in `plot` and `as_revman` methods for ReitsmaSubgroup models when 
  sensitivity and/or specificity are constrained to be common across subgroups.

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