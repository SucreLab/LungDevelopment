get_lower_ci <- function(data){
  b_ci_out <- suppressWarnings(boot.ci(boot(data = data,statistic = function(x,i) median(x[i]),R = 1000)))
  return(b_ci_out$percent[length(b_ci_out$percent) - 1])
}
get_upper_ci <- function(data){
  b_ci_out <- suppressWarnings(boot.ci(boot(data = data,statistic = function(x,i) median(x[i]),R = 1000)))
  return(b_ci_out$percent[length(b_ci_out$percent)])
}