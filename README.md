# PartMI <a href="https://danymukesha.github.io/PartMI/"><img src="man/figures/logo.png" align="right" align="right" height="90" alt="PartMI website" /></a>

{Part Mutual Information for Direct Association Networks} 

<!--
[![R-CMD-check](https://github.com/danymukesha/PartMI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danymukesha/PartMI/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/danymukesha/PartMI/branch/main/graph/badge.svg)](https://codecov.io/gh/danymukesha/PartMI)
-->

## Overview

**PartMI** implements Part Mutual Information (PMI) for 
inferring **direct associations** in biological networks. 
Unlike standard correlation or mutual information (MI), 
PMI conditions on background variables to remove **indirect dependencies**, 
substantially reducing false-positive edges.

### Why PMI?

| Method | Detects non-linear relationships | Removes indirect edges | Computationally feasible |
|---|---|---|---|
| Pearson correlation | No | No | Yes |
| Mutual Information (MI) | Yes | No | Yes |
| Partial correlation | No | Yes | Yes |
| **Part Mutual Information (PMI)** | **Yes** | **Yes** | **Yes** |

## Installation

```r
# install.packages("remotes")
remotes::install_github("danymukesha/PartMI")
```

## Quick Start

```r
library(PartMI)

# Simulate expression data
set.seed(42)
n <- 200
z <- rnorm(n)                    # hidden regulator
x <- 0.7 * z + rnorm(n, 0, 0.5) # gene X
y <- 0.7 * z + rnorm(n, 0, 0.5) # gene Y (indirect via Z)

# MI picks up the indirect association
mi(x, y)  # > 0 (false positive)

# PMI removes the common cause Z
pmi(x, y, z)  # ~ 0 (correct)
```

## Core Functions

| Function | Description |
|---|---|
| `mi()` | Mutual information (binning or kNN estimator) |
| `pmi()` | Part/conditional mutual information X vs Y given Z |
| `pmi_network()` | Build full network with significance testing |
| `normalize_data()` | Z-score, rank, minmax, or quantile normalization |
| `discretize_data()` | Quantile or equal-width discretization |

## References

- Zhao et al. (2016) *Part mutual information for quantifying direct associations in networks.* PNAS 113(18):5130-5135.
- Frenzel & Pompe (2007) *Partial Mutual Information for Coupling Analysis of Multivariate Time Series.* Phys. Rev. Lett. 99:204101.
- Kraskov, Stögbauer & Grassberger (2004) *Estimating mutual information.* Phys. Rev. E 69:066138.
