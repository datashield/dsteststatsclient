\name{chisq.test.combine}
\alias{chisq.test.combine}
\title{Runs Pearson's Chi-squared Test on non-pooled data}
\usage{
  chisq.test.combine(opals, a, b)
}
\arguments{
  \item{a}{}

  \item{b}{}
}
\value{
  statistic: the value the chi-squared test statistic.

  parameter: the degrees of freedom of the approximate
  chi-squared distribution of the test statistic

  p.value: the p-value for the test.

  method: a character string indicating the type of test
  performed, and whether Monte Carlo simulation or
  continuity correction was used.

  data.name: a character string giving the name(s) of the
  data.

  observed: the observed counts.

  expected: the expected counts under the null hypothesis.

  residuals: the Pearson residuals

  stdres: standardized residuals
}
\description{
  Performs Pearson's Chi-squared Test for Count Data
}
\examples{
\dontrun{
# put example here
}
}
\author{
  Gaye, A.
}

