# dtametaTMB

**Diagnostic Test Accuracy Meta-Analysis using Template Model Builder (TMB)**

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

- Frequentist estimation using binomial likelihoods  
- Support for **multiple thresholds per study** (Hoyer model)  
- Unified interface across DTA meta-analysis models  
- Summary ROC (SROC / HSROC) plots    
- Coupled forest plots  

---

## Installation

```r
# Install from GitHub (development version)
install.packages("devtools")
devtools::install_github("yourusername/dtametaTMB")
```