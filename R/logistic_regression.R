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
