#' Logistic regression modeling in a RDS study
#'
#' @param data Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest
#' @param formula An object of class 'formula'; description of the model with dependent and independent variables. Note that for linear regression the dependent variable must be continuous. Use reg.est.lm to perform regression with continuous dependent variable.
#' @param include.degree Logical; when TRUE the function performs weighted logistic regression with DEGREE used as weights.
#'
#' @return A data frame consisting of the point estimates and associated standard errors
#' @examples
#' reg_glm_out = reg.est.glm(data, formula, include.degree=F)
#' @export


reg.est.glm = function(data, formula, include.degree = F){

  data$inverse_weight = ave(data$DEGREE, data$boot_n, FUN = function(x) 1/x)
  df_groups = split(data, data$boot_n)


  if (include.degree == F) {
    point.estimates = lapply(df_groups, function(x) {
      models = glm(formula, data=x, family=quasibinomial(link='logit'))
      point.est = summary(models)$coefficients[,c(1)]
      return(point.est)
    })

    standard.errors = lapply(df_groups, function(x) {
      models = glm(formula, data=x, family=quasibinomial(link='logit'))
      point.est = summary(models)$coefficients[,c(2)]
      return(point.est)
    })
  }

  else {

    point.estimates = lapply(df_groups, function(x) {

      sv.design =  survey::svydesign(ids = ~0, # no clusters
                                     strata = NULL, # no strata
                                     data = x,
                                     weights = ~inverse_weight)

      models = survey::svyglm(formula, data=x, family=quasibinomial(link='logit'),
                              weights = inverse_weight, design = sv.design)
      point.est = summary(models)$coefficients[,c(1)]
      return(point.est)
    })

    standard.errors = lapply(df_groups, function(x) {
      sv.design =  survey::svydesign(ids = ~0, # no clusters
                                     strata = NULL, # no strata
                                     data = x,
                                     weights = ~inverse_weight)
      models = survey::svyglm(formula, data=x, family=quasibinomial(link='logit'),
                              weights = inverse_weight, design = sv.design)
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
