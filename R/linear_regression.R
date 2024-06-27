#' Linear regression modeling in an RDS study
#'
#' @param data Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest
#' @param formula An object of class 'formula'; description of the model with dependent and independent variables. Note that for linear regression the dependent variable must be continuous. Use reg.est.glm to perform regression with binary categorical dependent variable.
#' @param include.degree Logical; when TRUE the function performs weighted linear regression with DEGREE used as weights
#'
#' @return A data frame consisting of the point estimates and associated standard errors
#'
#' @examples
#' res_lm_out = reg.est.lm(data, formula, include.degree=G)
#'
#' @export


reg.est.lm = function(data, formula, include.degree = F){

  data$inverse_weight <- ave(data$DEGREE, data$boot_n, FUN = function(x) 1/x)
  df_groups = split(data, data$boot_n)


  if (include.degree == F) {
    point.estimates = lapply(df_groups, function(x) {
      models = lm(formula, data=x)
      point.est = summary(models)$coefficients[,c(1)]
      return(point.est)
    })

    standard.errors = lapply(df_groups, function(x) {
      models = lm(formula, data=x)
      point.est = summary(models)$coefficients[,c(2)]
      return(point.est)
    })
  }

  else if (include.degree == T){

    point.estimates = lapply(df_groups, function(x) {
      models = lm(formula, data=x, weight=inverse_weight)
      point.est = summary(models)$coefficients[,c(1)]
      return(point.est)
    })

    standard.errors = lapply(df_groups, function(x) {
      models = lm(formula, data=x, weight=inverse_weight)
      point.est = summary(models)$coefficients[,c(2)]
      return(point.est)
    })
  }

  pe = t(as.data.frame(point.estimates))
  se = t(as.data.frame(standard.errors))
  colnames(se)[1:ncol(se)] = paste('se_',colnames(se), sep='')
  results = apply(cbind(pe, se), 2, mean)
  return(results)

}
