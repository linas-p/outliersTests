#' Calculation of BP statistic
#'
#' The calculation of BP (Bagdonavičius-Petkevičius) statistic
#' for testing the hypothesis of absence of outliers
#' @param data A given data
#' @param distribution The distribution name (default: normal distribution)
#' @param alternative The alternative (default: two.sided)
#' @param s The number of the most remote z-scores used in the first step of outlier search (default: s=5)
#' @param location Location parameter (if it is specified). If location parameter is not specified, then location=NULL (default) and a robust estimate is calculated.
#' @param scale Scale parameter (if it is specified). If scale parameter is not specified, then  scale=NULL (default) and a robust estimate is calculated.
#' @return The list containing main statistic U, statistics U_i, ID of outliers, estimates location and scale
#' @keywords BP statistic
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' bp_statistic(x)
#'
bp_statistic <-  function(data, distribution = "norm",
  alternative = "greater", s = 5,
  location = NULL, scale = NULL) {
  n <- length(data);
  if (!(distribution %in% c("norm", "logis", "cauchy", "gumbel",
    "laplace", "lnorm", "extII"))){
    stop(paste("given distribution is supported yet"));
  }



  if ((s < 1) || (s > round(n/2))) {
    stop(paste("Statistics limit s is not valid"));
  }

  if(!(alternative %in% c("greater", "less", "two.sided"))){
    stop(paste("Choose alternative from greater/less/two.sided"));
  }

  order <- order(data);

  if (is.null(location) || is.null(scale)) {
    est <- get_robust_estimates(data, distribution);
  } else {
    est <- c()
    est$location <- location
    est$scale <- scale
  }

  bn <- get_bn(n, est$location, est$scale, distribution = distribution, alternative = alternative);
  an <- get_an(n, est$location, est$scale, distribution = distribution, alternative = alternative);

  if (alternative == "two.sided") {
    if (distribution == "gumbel" || distribution == "extII") {
      order <- order(data, decreasing = TRUE);
      idx1 <- order[1:s];
      idx2 <- order[n:(n-s+1)];
      rk1 <- ((data[idx1]) - bn[1])/an[1];
      rk2 <- (bn[2] + (data[idx2]))/an[2];
    } else {
      y <- abs(data)
      order <- order(y, decreasing = TRUE);
      idx <- order[1:s];
      rk <- (abs(y[idx]) - bn)/an;
    }
  } else if(alternative == "greater") {
    order <- order(data, decreasing = TRUE);
    idx <- order[1:s];
    rk <- (data[idx] - bn)/an;
  } else if(alternative == "less") {
    order <- order(data);
    idx <- order[1:s];
    rk <- -(data[idx] + bn)/an;
  }


  if(distribution == "cauchy") {
    rk <- 1-pchisq(2/(1+rk), 2*(1:s));
  } else if((distribution == "gumbel" | distribution == "extII")& alternative == "two.sided") {
      rk1 <- 1-pchisq(2*exp(-rk1), 2*(1:s));
      rk2 <- 1-pchisq(2*exp(rk2), 2*(1:s));
  } else {
    rk <- 1-pchisq(2*exp(-rk), 2*(1:s));
  }

  if((distribution == "gumbel" | distribution == "extII") & alternative == "two.sided") {
    return(list(Test_statistic_U = max(rk1), U_i = rk1, order = idx1, Ms2 = max(rk2), U_i2 = rk2, order2 = idx2,
                location = est$location, scale = est$scale))
  } else {
    return(list(Test_statistic_U = max(rk), U_i = rk, order = idx, location = est$location, scale = est$scale))
  }
}




#' BP test for outliers
#'
#' The calculation of the test statistic U and the p-value of the BP test for
#' outliers and finding of observations declared by the test as outliers.
#' @param data A given data
#' @param alpha The significance level
#' @param distribution The distribution name (default: normal distribution)
#' @param alternative The alternatives "two.sided"/"greater"/"less" (default: two.sided)
#' @param pvalue The indicator to return p-value (TRUE or FALSE (default)).
#'               The p-value is calculated using simulations SO IT MIGHT TAKE SOME TIME
#' @param pvalue_gen The number of generations for pvalue calculation; default=10000 takes up to 5 minutes
#' @return The list containing main statistic U, statistics U_i,  ID of outliers, estimates loc and scale,
#' and number of identified outliers
#' @keywords BP test
#' @export
#' @examples
#' set.seed(12)
#' x <- rnorm(100)
#' ks.test(x, "pnorm") # check normality
#' x[12:20] <- 4
#' ks.test(x, "pnorm") # check normality, not normal data
#' bp <- bp_test(x, pvalue = TRUE)
#' bp
#' x_after <- x[!bp$outlier]
#' ks.test(x_after, "pnorm") # check normality, after outliers removal data normal again
#'
#' set.seed(12)
#' x <- rcauchy(100)
#' ks.test(x, "pcauchy") # check cauchy distribution
#' x[12:22] <- 500
#' ks.test(x, "pcauchy") # check cauchy, not cauchy data
#' bp <- bp_test(x, distribution = "cauchy", pvalue = TRUE)
#' bp
#' x_after <- x[!bp$outlier]
#' ks.test(x_after, "pcauchy") # check cauchy, after outliers removal data cauchy again
#'
bp_test <-  function(data, alpha = 0.05, distribution = "norm",
  alternative = "greater", pvalue = FALSE, pvalue_gen = 1e4) {

  n <- length(data);
  if (!(distribution %in% c("norm", "logis", "cauchy", "gumbel",
    "laplace", "lnorm", "extII"))){
    stop(paste("given distribution is supported yet"));
  }

  if (!(alternative %in% c("less", "greater"))){
    stop(paste("given alternative is supported yet"));
  }

  if(!(alternative %in% c("greater", "less", "two.sided"))){
    stop(paste("Choose alternative from greater/less/two.sided"));
  }


  found_outliers <- FALSE;
  id <- c();
  statistic <- bp_statistic(data, distribution = distribution, alternative = alternative)
  ms <- statistic$Test_statistic_U;
  idx <- 1:n;

  outlier <- c();
  s0 <- 5

  critical <- get_critical(n, 5, alpha, distribution, alternative);

  if(distribution == "gumbel" | alternative == "two.sided") {
    found_outliers <- any((statistic$U_i > critical) | (statistic$U_i2 > critical));

    if(found_outliers == TRUE) {
      while( ((tail(statistic$U_i > critical, 1) == TRUE) | (tail(statistic$U_i2 > critical, 1) == TRUE))  & (length(data)*2 > n)) {
        statistic <- bp_statistic(data, distribution = distribution, alternative = alternative)
        critical <- get_critical(n, 5, alpha, distribution, alternative);

        found_local <- any((statistic$U_i > critical) | (statistic$U_i2 > critical));
        if(found_local) {
          if(statistic$Test_statistic_U > statistic$Test_statistic_U2 & statistic$Test_statistic_U > critical) { # reject right
            outlier <- c(outlier, idx[statistic$order[1]]);
            data <- data[-c(statistic$order[1])];
            idx <- idx[-c(statistic$order[1])]
          } else if(statistic$Test_statistic_U2 > critical) { # reject right
            nk <- length(data)
            outlier <- c(outlier, idx[statistic$order2[1]]);
            data <- data[-c(statistic$order2[1])];
            idx <- idx[-c(statistic$order2[1])]
          }
          #s <- tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)-1;
        }
      }

      statistic <- bp_statistic(data, distribution = distribution, alternative = alternative)
      critical <- get_critical(n, 5, alpha, distribution, alternative);
      if(length(which(as.numeric(statistic$U_i >= critical) != 0)) > 0) {
        outlier <- c(outlier, idx[statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)]]);
        data <- data[-c(statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)])];
        idx <- idx[-statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)]]
      }

      statistic <- bp_statistic(data, distribution = distribution, alternative = alternative)
      critical <- get_critical(n, 5, alpha, distribution, alternative);
      if(length(which(as.numeric(statistic$U_i2 >= critical) != 0)) > 0) {
        outlier <- c(outlier, idx[statistic$order2[1:tail(which(as.numeric(statistic$U_i2 >= critical) != 0), 1)]]);
        data <- data[-c(statistic$order2[1:tail(which(as.numeric(statistic$U_i2 >= critical) != 0), 1)])];
        idx <- idx[-statistic$order2[1:tail(which(as.numeric(statistic$U_i2 >= critical) != 0), 1)]]
      }


    }

  } else {
    found_outliers <- any((statistic$U_i > critical));
    if(found_outliers == TRUE) {
      while(tail(statistic$U_i > critical, 1) == TRUE & (length(data)*2 > n)) {
        statistic <- bp_statistic(data, distribution = distribution, alternative = alternative, location = statistic$location, scale = statistic$scale)
        n0 <- length(data)
        critical <- get_critical(n0, 5, alpha, distribution, alternative);
        found_local <- any((statistic$U_i > critical));
        if(found_local) {
          outlier <- c(outlier, idx[statistic$order[1]]);
          data <- data[-c(statistic$order[1])];
          idx <- idx[-c(statistic$order[1])]
        }
      }

      statistic <- bp_statistic(data, distribution = distribution, alternative = alternative, location = statistic$location, scale = statistic$scale)
      n0 <- length(data)
      critical <- get_critical(n0, 5, alpha, distribution, alternative);
      if(length(which(as.numeric(statistic$U_i >= critical) != 0)) > 0) {
        outlier <- c(outlier, idx[statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)]]);
        data <- data[-c(statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)])];
        idx <- idx[-statistic$order[1:tail(which(as.numeric(statistic$U_i >= critical) != 0), 1)]]
      }
    }
  }



  which_idexes <- rep(FALSE, n);

  if(length(outlier) > 0){
    which_idexes[outlier] <- TRUE;
  }

  if(pvalue) {
    p.value <- get_pvalue(ms, n, distribution = distribution, alternative = alternative, gen_num = pvalue_gen)
    return(list(found_outliers = found_outliers, outlier = which_idexes, id = outlier, p.value = p.value, number_of_outliers = length(outlier)));
  } else {
    return(list(found_outliers = found_outliers, outlier = which_idexes, id = outlier, number_of_outliers = length(outlier)));
  }



}


BP2TestReg <-  function(sample, distribution = "norm", alternative = "two.sided", correction = FALSE, s = 0) {
    n <- length(sample);
    if(s == 0) {
        s <- 5;
    }
    order <- order(sample)

    bn <- get_bn(n, distribution = distribution, alternative=alternative);
    an <- get_an(n, distribution = distribution, alternative=alternative);

    if (alternative == "two.sided") {
        y <- abs(sample)
        order <- order(y, decreasing = TRUE);
        idx <- order[1:s];
        rk <- (abs(y[idx]) - bn)/an;
    } else if(alternative == "greater") {
        order <- order(sample, decreasing = TRUE);
        idx <- order[1:s];
        rk <- (sample[idx] - bn)/an;
    } else if(alternative == "less") {
        order <- order(sample);
        idx <- order[1:s];
        rk <- -(sample[idx] + bn)/an;
    }


    rk <- 1-pchisq(2*exp(-rk), 2*(1:s));
    if(correction) {
        for(m in 1:s) {
            pars <- unifs[unifs[,1] == n & unifs[,2] == m, 3:4]
            rk[m] <- pbeta(rk[m], pars[1], pars[2]);
        }
    }

    return(list(Ms = max(rk), Ts = rk, order = idx));
}



BP2Test <-  function(sample, distribution = "norm", alternative = "greater", correction = FALSE,
                     s = 5) {
  n <- length(sample);
  order <- order(sample)

  est <- get_robust_estimates(sample, distribution);

  bn <- get_bn(n, est$location, est$scale, distribution = distribution, alternative = alternative);
  an <- get_an(n, est$location, est$scale, distribution = distribution, alternative = alternative);


  if (alternative == "two.sided") {
    y <- abs(sample)
    order <- order(y, decreasing = TRUE);
    idx <- order[1:s];
    rk <- (abs(y[idx]) - bn)/an;
  } else if(alternative == "greater") {
    order <- order(sample, decreasing = TRUE);
    idx <- order[1:s];
    rk <- (sample[idx] - bn)/an;
  } else if(alternative == "less") {
    order <- order(sample);
    idx <- order[1:s];
    rk <- -(sample[idx] + bn)/an;
  }
  rk <- 1-pchisq(2*exp(-rk), 2*(1:s));

  return(list(Ms = max(rk), Ts = rk, order = idx, mu = est$location, s = est$scale))
}



BP_regression_test <- function(sample, critical = NULL, distribution = "norm", n = length(sample), alternative = "two.sided",
                    alpha = 0.05, correction = FALSE, riba = 5, s = 5) {

    sample0 <- sample
    found <- FALSE;
    id <- c();
    statistic <- BP2Test(sample, distribution = distribution,
                         alternative = alternative, correction = correction,
                         s = 5)

    idx <- 1:n;

    outlier <- c();
    s0 <- s

    critical <- get_critical(s);
    found <- any((statistic$Ts > critical));

    if(found) {
        while(tail(statistic$Ts > critical, 1) == TRUE) {
            s <- s + 1;
            critical <- get_critical(s);
            statistic <- BP2Test(sample, distribution = distribution,
                                 alternative = alternative, correction = correction,
                                 s = s)
        }

        found_local <- TRUE;
        while(s > (riba-1) & found_local) { # not finished all s is "suspected"
            statistic <- BP2Test(sample, distribution = distribution,
                                 alternative = alternative, correction = correction,
                                 s = s)
            critical <- get_critical(s);
            found_local <- any((statistic$Ts > critical));
            if(found_local) {
                outlier <- c(outlier, idx[statistic$order[1]]);
                sample <- sample[-c(statistic$order[1])];
                idx <- idx[-c(statistic$order[1])]

                s <- tail(which(as.numeric(statistic$Ts >= critical) != 0), 1)-1;
            }
        }

        statistic <- BP2Test(sample, distribution = distribution,
                             alternative = alternative, correction = correction,
                             s = s)
        critical <- get_critical(s);
        if(length(which(as.numeric(statistic$Ts >= critical) != 0)) > 0) {
            outlier <- c(outlier, idx[statistic$order[1:tail(which(as.numeric(statistic$Ts >= critical) != 0), 1)]]);
            sample <- sample[-c(statistic$order[1:tail(which(as.numeric(statistic$Ts >= critical) != 0), 1)])];
            idx <- idx[-statistic$order[1:tail(which(as.numeric(statistic$Ts >= critical) != 0), 1)]]
        }
    }
    p <- rep(FALSE, n);

    if(length(outlier) > 0){
        p[outlier] <- TRUE;
    }

    return(list(found = found, outlier = p, id = outlier));

}


get_betas <- function(x, y, mode = "LTS", alternative = "two.sided", prop = NULL) {#, correction = FALSE

  if(mode == "OLS") {
    m <- lm(y ~ x);
  } else if(mode == "LTS") {
    dat <- as.data.frame(cbind(x, y))
    if (is.null(prop)) {
      prop <- 0.5
    }
    m <- ltsReg(y ~ .,  data = dat, alpha = prop); #alpha = get_k_n(length(y)),
  } else if(mode == "robust") {
    m <- rlm(y ~ x);
  } else if(mode == "LMSTS") {
    m <- ltsreg(y ~ x, method = "lts");
  } else if(mode == "LMSQS") {
    m <- lqs(y ~ x, method = "lqs");
  } else if(mode == "LMSMS") {
    m <- lmsreg(y ~ x, method = "lms");
  } else if(mode == "QR") {
    m <- rq(y ~ x);
  } else if(mode == "ROB") {
    m <- lmrob(y ~ x);
  }

  if(mode == "OLS") {
    hatsigma <- get_robust_estimates(m$residuals, type = "ML", alternative = alternative)$sigma;
    hatmu <- 0;
  } else {
    hatsigma <- m$scale;#est$sigma;
    hatmu <- 0; #est$mu;
  }
  #
  #
  if(is.null(dim(x))){
    n <- length(x);
  }else {
    n <- dim(x)[1];
  }

  ones <- rep(1, n);
  Z <- m$X;#cbind(ones, x);
  H <- Z %*% solve(t(Z) %*% Z) %*% t(Z);
  hii <- diag(H);
  d <- sqrt(1 - hii);

  ri <- (m$residuals) / (d*hatsigma); #/cl


  return(list(sigma = hatsigma, mu = hatmu, params = m$coefficients, ri = ri, hi = hii, ei = (m$residuals)/hatsigma));
}
