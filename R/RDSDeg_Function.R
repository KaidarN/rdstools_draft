Est_Deg<-function(sample,pPAT,covariates,maxiter,beta_restrict,threshold){
  # Actual degree (ActDeg) initiation
  sample<-sample%>%rowwise()%>%
    mutate(ActDeg=max(rpois(1,RepDeg), YPat,1))
  alpha<-c(log(mean(sample$RepDeg)),rep(0,length(covariates)))
  
  # Reporting behavior (RepDegBV) initiation
  sample$RepDegBV<-sample(c("error","exact","heap"),nrow(sample),replace = T)
  sample$RepDegBV<-factor(sample$RepDegBV,levels = c("error","exact","heap"))
  fit_beta<-multinom(data=sample,RepDegBV~ActDeg)
  beta<-coef(fit_beta)%>%as.matrix()
  beta<-rbind(beta,rep(0,2))
  rownames(beta)[3]<-"error"
  
  # Parameters of reporting rules initiation
  HeapUp<-sample%>%filter(RepDegBV=="heap")%>%
    mutate(grp=ifelse(RecRate<mean(sample$RecRate),"< ave. rate",">= ave. rate"))%>%
    mutate(k=(RepDeg-floor(ActDeg/5)*5)/5)%>%
    filter(k>0)
  freqUp<-HeapUp%>%group_by(grp,k)%>%count()
  
  delta_1<-exp(coef(lm(data=freqUp[freqUp$grp=="< ave. rate",],log(n)~k,weights = n))[2])
  
  delta_2_2<-exp(coef(lm(data=freqUp[freqUp$grp==">= ave. rate",],weights = n,log(n)~abs(k)))[2])
  
  HeapDown<-sample%>%filter(RepDegBV=="heap")%>%
    filter(RecRate>=mean(sample$RecRate))%>%
    mutate(k=(RepDeg-floor(ActDeg/5)*5)/5)%>%
    filter(k<=0)%>%
    mutate(K=floor(ActDeg/5)-1)%>%
    filter(K>0)
  wt<-HeapDown%>%group_by(K)%>%count()
  freqDown<-HeapDown%>%group_by(K,k)%>%count()
  rm_K<-unique(freqDown$K)[sapply(unique(freqDown$K),function(x)sum(freqDown$K==x)==1)]
  freqDown<-freqDown[!freqDown$K%in%rm_K,]
  wt<-wt[!wt$K%in%rm_K,]
  delta_2_1<-0
  for(i in unique(freqDown$K)){
    q<-exp(coef(lm(data=freqDown[freqDown$K==i,],weights = n,log(n)~abs(k)))[2])
    if(q>=1)
      q<-.99
    if(q<=0)
      q<-.01
    delta_2_1<-delta_2_1+wt$n[match(i,wt$K)]/sum(wt$n)*(1-q)
  }
  
  if(is.na(delta_2_1)){
    delta_2_1<-0
    w_1<-0
  }
  else
    w_1<-sum(wt$n)
  if(is.na(delta_2_2)){
    delta_2_2<-0
    w_2<-0
  }
  else
    w_2<-sum(freqUp[freqUp$grp==">= ave. rate","n"])
  
  delta_2<-delta_2_1*w_1/(w_1+w_2)+ delta_2_2*w_2/(w_1+w_2)
  
  if(delta_1<=0)
    delta_1<-.01
  if(delta_1>=1)
    delta_1<-.99
  if(delta_2<=0)
    delta_2<-.01
  if(delta_2>=1)
    delta_2<-.99
  X<-sample%>%
    mutate(Intercept=1)%>%
    dplyr::select(Intercept,covariates)%>%
    mutate_if(is.factor, ~ as.numeric(as.character(.x)))%>%
    as.matrix()
  out<-runEM(1:max(sample$RepDeg),maxiter,log(beta_restrict),threshold,
             alpha,beta,pPAT,X,sample$YPat,sample$RepDeg,
             sample$RecRate,mean(sample$RecRate[sample$RecRate>=0]), delta_1,delta_2,1)
  return(out)
  
}
