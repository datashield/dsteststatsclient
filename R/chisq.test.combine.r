#' Performs Pearson's Chi-squared Test for Count Data
#'
#' @title Runs Pearson's Chi-squared Test on non-pooled data
#'
#' @param a
#' @param b
#'
#' @return statistic: the value the chi-squared test statistic.
#' @return parameter: the degrees of freedom of the approximate chi-squared
#'            distribution of the test statistic
#' @return p.value: the p-value for the test.
#' @return method: a character string indicating the type of test performed, and
#'            whether Monte Carlo simulation or continuity correction was used.
#' @return data.name: a character string giving the name(s) of the data.
#' @return observed: the observed counts.
#' @return expected: the expected counts under the null hypothesis.
#' @return residuals: the Pearson residuals
#' @return stdres: standardized residuals
#' @author Gaye, A.
#' @examples
#' \dontrun{
#' # put example here
#'}
#' @export
#'
chisq.test.combine <- function (opals, a, b){
  # use the function table.2d.combine to create a 2d table of counts
  t2d <- datashield.table.2d.combine(opals,a,b)
  # use R function chisq.test to run the test
  results <- chisq.test(t2d)
  return(results)
}