# outliersTests

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/linas-p/outliersTests.svg?branch=master)](https://travis-ci.org/linas-p/outliersTests)

<!-- badges: end -->

The R package for test proposed in article [Multiple outlier detection tests for parametric models](https://www.mdpi.com/2227-7390/8/12/2156).

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

Simple method to check outliers for multiple distributions `investigate_sample`:

```R
investigate_sample(example1)
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

## Cite:

    Bagdonavičius, V.; Petkevičius, L. Multiple Outlier Detection Tests for Parametric Models. Mathematics 2020, 8, 2156.

or
        
     @article{bagdonaviciusmulti2020,
        title={Multiple Outlier Detection Tests for Parametric Models},
        author={Bagdonavičius, Vilijandas and Petkevičius, Linas},
        journal={Mathematics},
        volume={8},
        number={12},
        pages={2156},
        year={2020},
        publisher={Multidisciplinary Digital Publishing Institute}
      }
