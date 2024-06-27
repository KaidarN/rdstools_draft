#' Correlation analysis in a RDS study
#'
#' @param x Numeric vector; a continuous variable x
#' @param y Numeric vector; a continuous variable y
#' @param data Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest
#'
#' @return A numeric scalar value indicating correlation coefficient between x and y
#'
#' @examples
#' cor_coef = cor.est(x, y, data)
#' @export

cor.est = function(x, y, data) {
  df_groups = split(data, data$boot_n)
  cor.boot = lapply(df_groups, function(z) {
    results = cor(x,y)})
  cor.boot = mean(t(as.data.frame(cor.boot)))
}
