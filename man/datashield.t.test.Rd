\name{datashield.t.test}
\alias{datashield.t.test}
\title{Runs a student's t-test on non-pooled data}
\usage{
  datashield.t.test(opals, x, y = NULL, alternative = 1,
    mu = 0, paired = FALSE, var.equal = FALSE,
    conf.level = 0.95)
}
\arguments{
  \item{opals}{character strings that represent the URL of
  the servers where the study datasets are stored.}

  \item{x}{a (non-empty) numeric vector of data values}

  \item{y}{an optional (non-empty) numeric vector of data
  values.}

  \item{alternative}{character string specifying the
  alternative hypothesis, must be one of "two.sided"
  (default), "greater" or "less". You can specify just the
  initial letter.}

  \item{mu}{a number indicating the true value of the mean
  (or difference in means if you are performing a two
  sample test).}

  \item{paired}{a logical indicating whether you want a
  paired t-test.}

  \item{var.equal}{a logical variable indicating whether to
  treat the two variances as being equal. If TRUE then the
  pooled variance is used to estimate the variance
  otherwise the Welch. (or Satterthwaite) approximation to
  the degrees of freedom is used.}

  \item{confidence}{level of the interval}
}
\value{
  conf.level coefficients estimates of the parameters of
  the fitted model

  a list containing the following elements:

  the value of the t-statistic

  parameter the degrees of freedom for the t-statistic

  p.value the p-value for the test

  conf.int a confidence interval for the mean appropriate
  to the specified alternative hypothesis

  estimate the estimated mean or difference in means
  depending on whether it was a one-sample test or a
  two-sample test

  null.value the specified hypothesized value of the mean
  or mean difference depending on whether it was a
  one-sample test or a two-sample test

  alternative a character string describing the alternative
  hypothesis

  method a character string indicating what type of t-test
  was performed
}
\description{
  Performs one and two sample t-tests on vectors of data.
}
\details{
  Summary statistics are obtained from each of the data
  sets that are located on the distinct computers/servers.
  And then grand means and variances are calculated. Those
  are used for performing t-test.
}
\examples{
\dontrun{
# put example here
}
}
\author{
  Isaeva J.
}

