#' Calculation of significance level alpha_n
#'
#' For a given sample size n and the significance level alpha (alpha is the
#' probability for minimum one observation to fall into the outlier region)
#' the probability alpha_n for an  individual observation  to fall into the
#' outlier region is calculated.
#' @param n A sample size
#' @param alpha A significance level, default value 0.05
#' @return alpha_n
#' @keywords alpha_n
#' @export
#' @examples
#' get_alpha_n(1) # returns alpha_n = 0.05 since sample size n = 1 and default value of alpha is 0.05
#'
#' get_alpha_n(1, 0.1) # returns alpha_n = 0.1 since sample size n = 1

get_alpha_n <- function(n, alpha = 0.05) {
  if(alpha <= 0 | 1 <= alpha){
    stop(paste("alpha not beteen 0 and 1"));
  }
  if(any(n <= 0)){
    stop(paste("Sample size n not positive"));
  }
  return(1 - (1 - alpha)^(1 / n));
}

#' Robust estimates for location-scale families of distributions
#'
#' The robust estimates of the scale (Q_n) and location parameters for
#' location-scale families of distributions.
#' Implemented distributions: normal, gumbel, cauchy, laplace, logistic
#' @param x A sample
#' @param distribution The distribution name
#' @return ret The list structure with location and scale -  the estimates of location and scale parameters.
#' @keywords robust estimator
#' @export
#' @examples
#' x <- rnorm(100)
#' estimates <- get_robust_estimates(x)
#' estimates$location
#' estimates$scale
#'
#' x <- rcauchy(100)
#' estimates <- get_robust_estimates(x, "cauchy")
#' estimates
#'
#' x <- rlogis(100)
#' estimates <- get_robust_estimates(x, "logis")
#' estimates
#'

get_robust_estimates <- function(x, distribution = "norm") {
  require("robustbase");

  if(!(distribution %in% c("norm", "logis", "cauchy", "gumbel",
    "laplace", "lnorm"))){
    stop(paste("given distribution is supported yet"));
  }

  n <- length(x);
  mu <- NA;
  sigma <- NA;

  if(distribution == "lnorm") {
    x <- log(x);
  }

  k_const <- switch (distribution,
    norm = 2.22191445,
    logis = 1.307883,
    cauchy = 1.2071068,
    gumbel = 1.957615,
    laplace = 1.9305030,
    lnorm = 2.22191445
  )

  sigma <- Qn(x, constant = k_const);
  med <- median(x);
  mu <- switch (distribution,
    norm = med,
    logis = med,
    cauchy = med,
    laplace = med,
    gumbel = med - 0.3665129*sigma,
    lnorm = med
  )

  ret <- list(location=mu, scale=sigma)
  return(ret);
}


#' constant b_n
#'
#' The calculation of the constant b_n necessary for BP test
#' @param n A sample size
#' @param m An estimate of location parameter mu (default value (theoretical) 0)
#' @param s An estimate of shape parameter sigma (default value (theoretical) 1)
#' @param distribution The distribution name
#' @param alternative The choice of alternative
#' @return ret The b_n value
#' @keywords b_n constant
#' @export
#' @examples
#' get_bn(100, distribution = "norm")
#'
get_bn <- function(n, m = 0, s = 1, distribution = "norm", alternative = "greater") {

  bn <- NA;
  if(!(alternative %in% c("greater", "less", "two.sided"))){
    stop(paste("Choose alternative from greater/less/two.sided"));
  }

  if(!(distribution %in% c("norm", "logis", "cauchy", "gumbel",
    "laplace", "lnorm"))){
    stop(paste("given distribution is supported yet"));
  }

  if(n < 3){
    stop(paste("Sample size is not sufficient"));
  }

  if(alternative == "greater" | alternative == "less") {
    if(distribution == "norm") {
      bn <- m + s*qnorm((n-1)/n);
    } else if(distribution == "gumbel") {
      bn <- m - s*log(-log(1-1/n));
    } else if(distribution == "logis") {
      bn <- m + s*log(n-1);
    } else if(distribution == "cauchy") {
      bn <- m + s*n/pi;
    } else if(distribution == "lnorm") {
      bn <- exp(m + s*qnorm((n-1)/n))
    } else if(distribution == "laplace") {
      bn <- m + s*log(n/2);
    }
  } else { ## two.sided alternative
    if(distribution == "norm") {
      bn <- m + s*qnorm(1-1/(2*n));
    } else if(distribution == "logis") {
      bn <- m + s*log(2*n-1);
    } else if(distribution == "cauchy") {
      bn <- m + s*2*n/pi;
    } else if(distribution == "laplace") {
      bn <- m + s*log(n);
    } else if(distribution == "lnorm") {
      bn <- m + s*log(n);
    } else if(distribution == "gumbel") {
      bn1 <- m + s*log(log(n));
      bn2 <- m - s*log(-log(1-1/n));
      bn <- c(bn2, bn1)
    }
  }
  return(bn);
}


#' constant a_n
#'
#' The calculation of the constant a_n nessisary for BP test
#' @param n A sample size
#' @param m An estimate of location parameter mu (default value (theoretical) 0)
#' @param s An estimate of shape parameter sigma (default value (theoretical) 1)
#' @param distribution The distribution name
#' @param alternative The choice of alternative
#' @return ret The a_n value
#' @keywords a_n constant
#' @export
#' @examples
#' get_an(100, distribution = "norm")
#'
#'
get_an <- function(n, m = 0, s=1, distribution = "norm", alternative = "greater"){
  an <- NA;

  if(!(alternative %in% c("greater", "less", "two.sided"))){
    stop(paste("Choose alternative from greater/less/two.sided"));
  }

  if(!(distribution %in% c("norm", "logis", "cauchy", "gumbel",
    "laplace", "lnorm"))){
    stop(paste("given distribution is supported yet"));
  }

  if(n < 3){
    stop(paste("Sample size is not sufficient"));
  }


  if(alternative == "greater" | alternative == "less") {
    if(distribution == "norm") {
      an <- s/qnorm(1 - 1/n)
    } else if(distribution == "gumbel") {
      an <- -s/((n-1)*log(1-1/n))
    } else if(distribution == "logis") {
      an <- s*n/(n-1);
    } else if(distribution == "cauchy") {
      an <- s*n/pi;
    } else if(distribution == "lnorm") {
      bn <- get_bn(n, m, s, distribution);
      an <- s*bn/(n*dnorm(qnorm(1 - 1/n)))
    } else if(distribution == "laplace") {
      an <- s;
    }
  } else {
    if(distribution == "norm") {
      an <- s/qnorm(1 - 1/(2*n))
    }  else if(distribution == "logis") {
      an <- s*4*n/(2*n-1);
    } else if(distribution == "cauchy") {
      an <- s*2*pi/(2*n*(sin(pi/(2*n)))^2);
    } else if(distribution == "laplace") {
      an <- 2*s;
    }  else if(distribution == "lnorm") {
      bn <- get_bn(n, m, s, distribution);
      an <- s*bn/(n*dnorm(qnorm(1 - 1/(2*n))))
    } else if(distribution == "gumbel") {
      an1 <- s/(log(n))
      an2 <- -s/((n-1)*log(1-1/n))
      an <- c(an2, an1)
    }
  }
  return(an);
}

#' Critical values of BP test
#'
#' The method with estimated aproximations of critical values.
#' The unimplemented cases using rough approximation instead.
#' @param n A sample size
#' @param s A upper limit (default s=5)
#' @param alpha The significance level  (default: 0.05)
#' @param distribution The distribution name
#' @param alternative The choice of alternative ("two.sided"/"greater"/"less")
#' @return ret The critical value of statistics U
#' @keywords BP critical values
#' @export
#' @examples
#' get_critical(100)
#'
#'
get_critical <- function(n, s = 5, alpha = 0.05, distribution = "norm", alter = "two.sided") {
  critical <- NA;

  if (alter == "two.sided") {
    if (distribution == "norm") {
      if(alpha == 0.05) {
        if (s == 5) {
          x <- log(log(log(log(n))));
          if(n %% 2 == 1) {
            critical <- 1 - (0.015929         -0.006473 * x)
          }else {
            critical <- 1 - (0.003551         -0.034848  * x)
          }
          return(critical)
        }
      }
    } else if (distribution == "logis") {
      if(alpha == 0.05) {
        if (s == 5) {
          x <- log(log(log(log(n))));
          if(n %% 2 == 1) {
            critical <- 1 - (0.010536      -0.006989  * x)
          }else {
            critical <- 1 - (0.005306      -0.018044  * x)
          }
          return(critical)
        }
      }
    } else if (distribution == "laplace") {
      if(alpha == 0.05) {
        if (s == 5) {
          x <- log(log(log(log(n))));
          if(n %% 2 == 1) {
            critical <- 1 - (0.009657      -0.006613 * x)
          }else {
            critical <- 1 - (0.006819      -0.013780 * x)
          }
          return(critical)
        }
      }
    } else if (distribution == "gumbel") {
      if(alpha == 0.05) {
        if (s == 5) {
          x <- log(log(log(log(n))));
          if(n %% 2 == 1) {
            critical <- 1 - (0.005755      -0.002166   * x)
          }else {
            critical <- 1 - ( 0.00302       -0.00840   * x)
          }
          return(critical)
        }
      }
    } else if (distribution == "cauchy") {
      if(alpha == 0.05) {
        if (s == 5) {
          x <- log(log(log(log(n))));
          if(n %% 2 == 1) {
            critical <- 1 - ( 0.009183      -0.003114  * x)
          }else {
            critical <- 1 - (0.008783      -0.004095   * x)
          }
          return(critical)
        }
      }
    }
  }

  if(is.na(critical)){
    warning("At given s, alpha and alternative, exact critical values is not implemented.
      Rough approximation is used insted.")
    critical <- (1-alpha)^(1/s)
  }
  return(critical)

}

#' Main BP statistic
#'
#' The calculation of the value of the BP test statistic U.
#' @param data A given data
#' @param distribution The distribution name
#' @param alternative The choice of alternative
#' @return ret The value of the statistic U
#' @keywords BP statistic
#' @export
#' @examples
#' x <- rnorm(100)
#' get_main_bp_statistic(x)
#'
#'
get_main_bp_statistic <-  function(data, distribution = "norm",
  alternative = "two.sided") {
  bp <- bp_statistic(data, distribution, alternative)
  return(bp$Test_statistic_U)
}

#' p-value of BP test
#'
#' The use simulation and for given conditions calculate exact p-value.
#' Since simulation is used it might take few minutes.
#' @param statistic_val A given BP statistics value
#' @param n The size of the sample
#' @param distribution The distribution name
#' @param alternative The chooise of alternative
#' @param gen_num The number of simulations to used in order to calculate p-value.
#' @return ret The p-value of given BP statistics
#' @keywords BP p-value
#' @export
#' @examples
#' x <- rnorm(100)
#' bp <- bp_statistic(x)
#' get_pvalue(bp$Test_statistic_U, 100)
#'
get_pvalue <- function(statistic_val, n, distribution = "norm", alternative = "two.sided", gen_num = 1e4) {
  r_samples <- switch(distribution,
    "norm" = lapply(rep(n, gen_num), rnorm),
    "cauchy" = lapply(rep(n, gen_num), rcauchy),
    "logis" = lapply(rep(n, gen_num), rlogis),
    "laplace" = lapply(rep(n, gen_num), rlaplace),
    "lnorm" = lapply(rep(n, gen_num), rlnorm),
    "gumbel" = lapply(rep(n, gen_num), rgumbel)
  )
  statistics <- lapply(r_samples, get_main_bp_statistic, distribution = distribution, alternative=alternative)
  p.value <- mean(unlist(statistics) > statistic_val)
  return(p.value)
}


#' outliersTests: A package containing statistical tests of identification
#' unknown number of outliers
#'
#' Statistical test for various location-scale family distributions
#' normal, logistic, cauchy, laplace to test does sample contain outliers and
#' identifying those outliers in sample.
#'
#' @section outliersTests functions:
#'
#'
#' \code{\link{bp_test}}: The calculation of the test statistic U and the p-value of the BP test for
#' outliers and finding of observations declared by the test as outliers.
#'
#' \code{\link{get_robust_estimates}}:  The robust estimates of the scale (Q_n) and location parameters for
#' location-scale families of distributions.
#' Implemented distributions: normal, gumbel, cauchy, laplace, logistic
#'
#' @section outliersTests package data:
#'
#' \code{\link{example1}}: - example data
#'
#'
#' @docType package
#' @name outliersTests
NULL
