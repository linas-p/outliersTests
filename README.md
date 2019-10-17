# outliersTests

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/linas-p/outliersTests.svg?branch=master)](https://travis-ci.org/linas-p/outliersTests)
[![Codecov test coverage](https://codecov.io/gh/linas-p/outliersTests/branch/master/graph/badge.svg)](https://codecov.io/gh/linas-p/outliersTests?branch=master)

<!-- badges: end -->



The purpose of outliersTests is to provide statistical tests for absence of outliers hypothesis for various location-scale families of  distributions and give procedures for outlier identification when the hypothesis is rejected.  Normal, logistic, Cauchy, Laplace, Gumbel families are included. In the case of shape-scale families such as Weibull, lognormal or  loglogistic it is sufficient to take logarithms of observations and apply all procedures for Gumbel, normal or logistic families, respectively. 

Key features:

* The most important functions: `bp_test`, `get_robust_estimates`.



## Installation

To get the current development version from github:

```R
# install.packages("devtools")
devtools::install_github("linas-p/outliersTests")
library(outliersTests)
```




## Usage

Most simple usage just call method `bp_test`:

```R
bp_test(example1)
```

For specified data more parameters can be indicated:

```R
set.seed(12)
x <- rcauchy(100)
ks.test(x, "pcauchy") # check cauchy distribution
x[12:22] <- 500
ks.test(x, "pcauchy") # check cauchy, not cauchy data
bp <- bp_test(x, distribution = "cauchy", pvalue = TRUE)
bp
x_after <- x[!bp$outlier]
ks.test(x_after, "pcauchy") # check cauchy, after outliers removal data cauchy again
```

The documentation of the usage is accesible as:

```R
?outliersTests # main description of package
?bp_test # documentation for the usage of the BP test
```
