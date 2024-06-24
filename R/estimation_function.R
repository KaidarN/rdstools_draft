rm(list=ls())
# hawk data
# toy data
library(tidyverse)
toy_data_rds = read_csv('C:/Users/kai_kasados/Desktop/RDS/package/toy_data_handle.csv')
toy_data_rds$age = rnorm(nrow(toy_data_rds), 40, 8)
## mean and median of bootstraps
## Hmisc has weighted variance
## matrixStats::weightedMedian() has weighted median

results$DEGREE = rnorm(nrow(results), 100, 10)

#descr.stats.est = function(data, x, include.degree = F){

#         data$inverse_weight <- ave(results$DEGREE, results$boot_n, FUN = function(x) 1/x) 
         
#if (include.degree == F){
#         mean.boot = mean(aggregate(x ~ data$boot_n, FUN = mean)[,2])
#         var.boot = mean(aggregate(x ~ data$boot_n, FUN = var)[,2])
#         median.boot = mean(aggregate(x ~ data$boot_n, FUN = median)[,2])
#         
#         ci.boot_lower = mean.boot-(2*sqrt(var.boot))
#         ci.boot_upper = mean.boot+(2*sqrt(var.boot))
#         results = as.data.frame(rbind(c('mean', 'median', 'ci_lower', 'ci_upper'),
#                                       c(round(mean.boot, 2), round(median.boot, 2), 
#                                         round(ci.boot_lower, 2), round(ci.boot_upper, 2))))
#         
#         colnames(results) <- results[1, ]
#         results = results[-1,]
#         
#}
         
#else {
#          df_groups = split(data, data$boot_n)  
#          
#         mean.boot = lapply(df_groups, function(x) {
#           res = weighted.mean(x, data$inverse_weight)
#           return(point.est)
#         })
#         mean.boot  = mean(aggregate(x~data$boot_n, data = data, FUN = w.mean(x, data$inverse_weight))[,2])
#         results = mean.boot
#}  
         
#         return(results)
#}

#res = descr.stats.est(results, results$data.cov, include.degree = T)
#res

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

res = descr.stats.est(results, results$data.cov)
res

results$x = rnorm(nrow(results), 10, 100)

reg.est.lm = function(data, formula, include.degree = F){
 
  data$inverse_weight <- ave(results$DEGREE, results$boot_n, FUN = function(x) 1/x) 
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

df = reg.est.lm(results, x~data.cov, include.degree = F)
df




results$gender =  sample(c(0, 1), size = nrow(results), replace=T)
res = glm(results$gender ~ results$x, data=results, family=quasibinomial(link='logit'))
summary(res)

### glm
###
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


cor.est = function(x, y, data) {
       df_groups = split(data, data$boot_n)
       cor.boot = lapply(df_groups, function(z) {
         results = cor(x,y)})
       cor.boot = mean(t(as.data.frame(cor.boot)))
}

res = cor.est(results$data.cov, results$data.cov, results)


res = t(as.data.frame(res))
mean(res)
res =  reg.est.glm(results, gender ~ x) 
res  

res = reg.est.glm(results, gender ~ x, include.degree = T)
res

