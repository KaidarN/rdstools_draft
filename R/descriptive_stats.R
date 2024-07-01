#' Descriptive statistics in an RDS study
#'
#' @param data Data frame; The output from bootstrap_RDS with boot_n and respondent ID
#' @param x Numeric vector; A variable of interest
#'
#' @return Data frame; Containes mean, median, variance and associated 95\% CIs
#'
#' @examples
#' data('RDStoydata')
#'
#' # Preprocess data with RDSdata function
#' rds_data <- RDSdata(data = RDStoydata,unique_id = "ID",
#' redeemed_coupon = "CouponR",
#' issued_coupon = c("Coupon1",
#'                  "Coupon2",
#'                  "Coupon3"),
#'                degree = "Degree",
#'                result = c('Age','Sex'))
#'
#'
#' # Run bootstrap_RDS with rds_data
#' results = bootstrap_RDS(RESPONDENT_ID = rds_data$ID, SEED_ID = rds_data$S_ID,
#' SEED = rds_data$SEED, RECRUITER_ID = rds_data$R_CP,
#' data.cov = rds_data[,c('Age', 'Sex')], type = 'boot_chain_one', n.times = 100,
#' return_data = T)
#'
#' # Calculate descriptive statistics using results from bootstrap_RDS
#' out = descr.stats.est(results, as.numeric(results$Age))
#' out
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
