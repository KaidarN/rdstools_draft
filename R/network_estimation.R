
#' Estimating individual degree size in RDS study.
#'
#' @param sample data.frame; An RDS sample, including information of:
#' ID (Respondent unique ID),
#' recruiterID (Recruiter's ID, not necessary for this function),
#' seedID (Seed ID, not necessary for this function),
#' wave (Wave number, not necessary for this function),
#' RepDeg (Respondent's reported degree),
#' YPat (Respondent's reported peers named Pat within target RDS population ),
#' RecRate (Respondent's recruitment rate, i.e., the ratio of the number of recruits to the number of issued coupons ),
#' Covariates (Subject-specific characteristics. If there exists continuous variables with large scales, we recommend a log- transform).
#' @param pPAT scalar; External knowledge of the  size proportion of subpopulation named "Pat" within the target RDS population.
#' @param covariates vector; A vector of subject-specific characteristics that are assumed to be associated with the true degree distribution via a Poisson regression model.
#' @param maxiter scalar; Maximum number of iterations used in EM algorithm. Default is 1000.
#' @param beta_restrict scalar; Bound of coefficients from the parametric regression model of the underlying reporting behavior. Default is 10000.
#' @param threshold scalar; Convergence threshold. Default is 0.005.
#' @return A list consisting of the following elements:
#' \item{degree}{vector; estimated individual true degree size.}
#'\item{alpha}{vector; Estimated coefficients obtained from the parametric regression model of the true degree
#'adjusted for subject-specific characteristics.}
#'\item{beta}{matrix; Estimated coefficients obtained from the parametric regression model of the underlying reporting behavior
#'adjusted for the individual true degree. }
#'\item{delta_1, delta_2, phi}{scalars; Estimated parameters obtained from models of reporting rules.}
#' @useDynLib RDStools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' out <- Est_Deg(sample,pPAT=0.01,
#'                covariates=c("Male","Black","Married","LowEduc","Age","Income","DETROIT",
#'                             "M_DQ3","M_DQ7","M_DQ29_SUPPORT","M_GQ11","M_HQ9","Duration"),
#'                maxiter=1000,beta_restrict=10000,threshold=0.005)
#'
#' @export


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
