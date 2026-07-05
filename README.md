# dtametaTMB

**Diagnostic Test Accuracy Meta-Analysis using Template Model Builder**

`dtametaTMB` provides a unified framework for frequentist meta-analysis of diagnostic
test accuracy (DTA) studies in R. It implements several widely used models
within a consistent interface, including:

- the Reitsma bivariate random-effects model  
- the Rutter and Gatsonis (HSROC) model  
- the Hoyer threshold-based bivariate time-to-event model  

The package uses **Template Model Builder (TMB)** for efficient likelihood-based
estimation of complex hierarchical models.

---

## Key Features

- Frequentist estimation using exact binomial likelihoods
- Reitsma bivariate random-effects model
- Rutter and Gatsonis (HSROC) model
- Hoyer multiple-threshold model
- Subgroup analyses and meta-regression
- Likelihood-ratio tests for nested model comparisons
- Support for parameter constraints in sparse-data settings
- Summary ROC (SROC/HSROC) plots
- Coupled forest plots
- Unified interface across DTA meta-analysis models
- Validation against examples from the Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy
