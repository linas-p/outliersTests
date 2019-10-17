# outliersTests


The aim of outliersTests is to provide a the statistical tests for various location-scale family distributions: normal, logistic, cauchy, laplace to test does sample contain outliers and identifying those outliers in sample.

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

If you have your on data you might used more parameters:

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
