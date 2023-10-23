# DEMIC

<!-- badges: start -->
  [![R-CMD-check](https://github.com/Ulthran/DEMIC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ulthran/DEMIC/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/Ulthran/DEMIC/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Ulthran/DEMIC?branch=master)
  [![DOI:10.1038/s41592-018-0182-0](https://badgen.net/badge/Published%20in/Nat%20Methods/blue)](https://doi.org/10.1038/s41592-018-0182-0)
<!-- badges: end -->

## Introduction



## Installation

The development version can be installed from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Ulthran/DEMIC")
```

## Basic Usage

```r
library(demic)
eptrs1 <- est_ptr(ContigCluster1)
eptrs2 <- est_ptr(ContigCluster2)
```

## Docs and Additional Help

This package is designed to be used with the COV3 output files of [PyCov3](https://github.com/Ulthran/pycov3). Its usage can be seen in the context of the sunbeam extension [sbx_demic](https://github.com/Ulthran/sbx_demic). Please cite Gao, Y., Li, H. Quantifying and comparing bacterial growth dynamics in multiple metagenomic samples. Nat Methods 15, 1041–1044 (2018). https://doi.org/10.1038/s41592-018-0182-0.
