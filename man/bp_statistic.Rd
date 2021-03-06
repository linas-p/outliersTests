% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BPTest.R
\name{bp_statistic}
\alias{bp_statistic}
\title{Calculation of BP statistic}
\usage{
bp_statistic(
  data,
  distribution = "norm",
  alternative = "two.sided",
  s = 5,
  location = NULL,
  scale = NULL
)
}
\arguments{
\item{data}{A given data}

\item{distribution}{The distribution name (default: normal distribution)}

\item{alternative}{The alternative (default: two.sided)}

\item{s}{The number of the most remote z-scores used in the first step of outlier search (default: s=5)}

\item{location}{Location parameter (if it is specified). If location parameter is not specified, then location=NULL (default) and a robust estimate is calculated.}

\item{scale}{Scale parameter (if it is specified). If scale parameter is not specified, then  scale=NULL (default) and a robust estimate is calculated.}
}
\value{
The list containing main statistic U, statistics U_i, ID of outliers, estimates location and scale
}
\description{
The calculation of BP (Bagdonavičius-Petkevičius) statistic
for testing the hypothesis of absence of outliers
}
\examples{
set.seed(123)
x <- rnorm(100)
bp_statistic(x)

}
\keyword{BP}
\keyword{statistic}
