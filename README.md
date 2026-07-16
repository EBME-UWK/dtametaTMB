# dtametaTMB

**Diagnostic Test Accuracy Meta-Analysis using Template Model Builder**

`dtametaTMB` provides a unified framework for frequentist meta-analysis of diagnostic
test accuracy (DTA) studies in R. It implements conventional gold-standard models,
latent class extensions for imperfect reference standards, subgroup analyses, and
multiple-threshold models within a consistent interface.

---

## Key Features

### Models

- Reitsma bivariate random-effects model
- Rutter and Gatsonis (HSROC) model
- Hoyer multiple-threshold model
- Latent class Reitsma and HSROC models for studies without a perfect reference standard

### Extensions

- Subgroup analyses
- HSROC meta-regression
- Parameter constraints for sparse-data settings
- Likelihood-ratio tests for nested model comparisons

### Output

- Summary ROC (SROC/HSROC) plots
- Coupled forest plots
- RevMan-compatible exports

### Implementation

- Frequentist estimation using exact binomial likelihoods
- Template Model Builder (TMB) backend
- Unified interface across model families
- Validation against published examples and Cochrane Handbook analyses