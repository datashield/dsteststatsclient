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
  # use the same approach as in the function table.2d.combine to create a 2d table of counts
  cally <- call("table.2d", a, b) 
  output.object <- datashield.aggregate(opals, cally)
  
  num.sources <- length(output.object)
  table.valid<-rep(NA,num.sources)
  
  for(j in 1:num.sources)
  {
    table.valid[j]<-(sum(output.object[[j]])>=0)
  }
  
  tab.final<-as.matrix(output.object[[1]])*0
  
  
  for(j in 1:num.sources)
  {
    if(table.valid[j]==1)
    {
      tab.final<-tab.final+as.matrix(output.object[[j]])
    }
  }
  
  
  total.count<-sum(tab.final)
  
  numrows<-dim(tab.final)[1]
  numcols<-dim(tab.final)[2]
  
  table.combination.valid<-1
  
  for(j in 1:numrows)
  {
    for(m in 2:num.sources)
    {
      if(dimnames(output.object[[m-1]])$Var1[j]!= dimnames(output.object[[m]])$Var1[j])
      {table.combination.valid<-0}
    }
  }
  
  for(j in 1:numcols)
  {
    for(m in 2:num.sources)
    {
      if(dimnames(output.object[[m-1]])$Var2[j]!= dimnames(output.object[[m]])$Var2[j])
      {table.combination.valid<-0}
    }
  }
  
  col.sum.vect<-rep(1,numrows)
  row.sum.vect<-rep(1,numcols)
  
  col.tots<-t(tab.final)%*%col.sum.vect
  row.tots<- tab.final%*%row.sum.vect
  
  tab.final.counts<-tab.final
  cat("\nThe 2D table used for chi-squared test\n\n")
  print(tab.final.counts)
  
  # use R function chisq.test to run the test
  results <- chisq.test(tab.final.counts)
  return(results)
}