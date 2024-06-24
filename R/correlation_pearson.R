cor.est = function(x, y, data) {
  df_groups = split(data, data$boot_n)
  cor.boot = lapply(df_groups, function(z) {
    results = cor(x,y)})
  cor.boot = mean(t(as.data.frame(cor.boot)))
}
