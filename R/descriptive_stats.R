#' Descriptive statistics in an RDS study
#'
#' @param data Data frame; The output from bootstrap_RDS with boot_n and respondent ID
#' @param x Numeric vector; A variable of interest
#'
#' @return Data frame; Containes mean, median, variance and associated 95\% CIs
#'
#' @examples
#' out = descr.stats.est(data, x)
#'
#' @export

descr.stats.est = function(data, x){

  mean.boot = mean(aggregate(x ~ data$boot_n, FUN = mean)[,2])
  var.boot = mean(aggregate(x ~ data$boot_n, FUN = var)[,2])
  median.boot = mean(aggregate(x ~ data$boot_n, FUN = median)[,2])

  ci.boot_lower = mean.boot-(2*sqrt(var.boot))
  ci.boot_upper = mean.boot+(2*sqrt(var.boot))
  results = as.data.frame(rbind(c('mean', 'median', 'ci_lower', 'ci_upper'),
                                c(round(mean.boot, 2), round(median.boot, 2),
                                  round(ci.boot_lower, 2), round(ci.boot_upper, 2))))

  colnames(results) <- results[1, ]
  results = results[-1,]
  return(results)
}
