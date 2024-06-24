#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/random.hpp>
#include <random>
#include <cmath>
#include <math.h>  
#include <RcppNumerical.h>
#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <limits>
using namespace Numer;
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo,RcppParallel,BH,RcppEigen,RcppNumerical)]]


class Mintegrand_func1: public Func
{
private:
  const double a;
  const double b;
public:
  Mintegrand_func1(double a_, double b_) : a(a_), b(b_) {}
  // calculate x^a*(1-x)^{b-1}*log(1-x)dx 
  double operator()(const double& x) const
  {
    return pow(x,a)*pow(1-x,b-1)*log(1-x) ;
  }
  
};
class Mintegrand_func2: public Func
{
private:
  const double a;
  const double b;
public:
  Mintegrand_func2(double a_, double b_) : a(a_), b(b_) {}
  // calculate  x^a*(1-x)^{b-1}*log(1-x)*log(1-x)dx 
  double operator()(const double& x) const
  {
    return pow(x,a)*pow(1-x,b-1)*log(1-x)*log(1-x) ;
  }
  
};
// [[Rcpp::export]]
double IntegralFunc1(double upper , double a, double b)
{  // calculate \int_{0}^{upper} x^a*(1-x)^{b-1}*log(1-x)dx 
  Mintegrand_func1 f(a, b);
  double err_est;
  int err_code;
  return integrate(f, 0.0, upper, err_est, err_code);
}
// [[Rcpp::export]]
double IntegralFunc2(double upper , double a, double b)
{
  // calculate \int_{0}^{upper} x^a*(1-x)^{b-1}*log(1-x)*log(1-x)dx 
  Mintegrand_func2 f(a, b);
  double err_est;
  int err_code;
  return integrate(f, 0.0, upper, err_est, err_code);
}
typedef std::vector<double> stdvec;
typedef std::vector< std::vector<double> > stdvecvec;
//[[Rcpp::export]]
stdvecvec mat_to_std_vec(mat &A) {
  //convert arma::mat to vec::ector
  stdvecvec V(A.n_rows);
  for (size_t i = 0; i < A.n_rows; ++i) {
    V[i] = conv_to< stdvec >::from(A.row(i));
  };
  return V;
}

//[[Rcpp::export]]
double logfactorial(int n) {
  if (n < 0) {
    // Factorial is not defined for negative numbers
    return std::numeric_limits<double>::quiet_NaN();
  } else if (n == 0) {
    // log(0!) = log(1) = 0
    return 0.0;
  } else {
    // For positive n, calculate the log factorial
    double result = 0.0;
    for (int i = 1; i <= n; ++i) {
      result += std::log(i);
    }
    return result;
  }
  
  // This should be unreachable, but it's provided as a safeguard.
  return std::numeric_limits<double>::infinity();
}
//[[Rcpp::export]]
vec cal_gamma(int di,mat beta){
  //calculate (gamma_exact,gamma_heap,gamma_guess)|(di,beta): 3 types of RepBV j, Pr(h_{ij}=1|di,beta) 
  colvec Ai={1,di};
  vec exp_elems (3);//output (gamma_exact,gamma_heap,gamma_guess)
  double x1=dot(Ai,beta.row(0));//beta_{exact,0}+beta_{exact,1}di
  double x2=dot(Ai,beta.row(1));//beta_{heap,0}+beta_{heap,1}di
  exp_elems(0)=x1;
  exp_elems(1)=x2;
  exp_elems(2)=0;
  exp_elems.for_each([](double& val){val=exp(val);});
  if(exp_elems.has_nan()==0&&exp_elems.has_inf()==0){
    double denom= accu(exp_elems);
    exp_elems=exp_elems/denom;
  }
  else{
    if((!isinf(exp(x2)))&&(!isnan(exp(x2)))){//x2 is valid
      if(x1>0)//exp(x1) is positive infinity
        exp_elems={1,0,0};
      else//exp(x1) is negative infinity
        exp_elems={0,exp(x2)/(1+exp(x2)),1/(1+exp(x2))};  
    }
    else{//x2 is NOT valid
      if((!isinf(exp(x1)))&&(!isnan(exp(x1)))){//x1 is valid
        if(x2>0)//exp(x2) is positive infinity
          exp_elems={0,1,0};
        else//exp(x2) is negative infinity
          exp_elems={exp(x1)/(1+exp(x1)),0,1/(1+exp(x1))};
      }
      else{//x1 and x2 both are NOT valid
        if(x1>0){
          if(x2>0)
            exp_elems={.5,0.5,0};
          else
            exp_elems={1,0,0};
        }
        else{
          if(x2<0)
            exp_elems={0,0,1};
          else
            exp_elems={0,1,0};
        }
      } 
    }
  }
  if(exp_elems.has_nan()|exp_elems.has_inf())
    Rcout<<"error in gamma";
  return exp_elems;
}

//[[Rcpp::export]]
arma::mat cal_complete( // here should be arma, why mat is not highlighted
    int di,
    double mu,
    arma::mat beta,
    const int yi,
    const int di_rep,
    const double recruit_rate,
    const double ave_rate,
    const double pPAT,
    const double delta_1,
    const double delta_2,
    const double phi
  ){

  //return individual-level Pr(y_obs,y_mis={di,h_{ij}=1})
  mat out(3,6); //(di,{h_{ij}}_{j=1,2,3},{gamma_j},Pr(y_obs,y_mis={di,h_{ij}=1}))
  out.zeros();
  
  out.col(0)=vec({di,di,di});out(0,1)=1;out(1,2)=1;out(2,3)=1;
  if(di<yi)//ActDeg must be >= number of PAT friends
    return out;
  
  double f_yi=exp(yi*log(pPAT)+(di-yi)*log(1-pPAT)+logfactorial(di)-logfactorial(yi)-logfactorial(di-yi));//yi~BIN(di,pPAT)
  double f_di=exp(di*log(mu)-mu-logfactorial(di))/(1-exp(-mu))*(di>0);//di~0-truncated Poisson(mu)
  vec gamma=cal_gamma(di,beta);//gamma_exact,gamma_heap,gamma_guess, which is Pr(h_{ij}=1)
  out.col(4)=gamma;
  //next calculate Pr(y_obs,y_mis={di,h_{ij}=1}))=f_yi*f_di*gamma*Pr(RepDeg|ActDeg,h_{ij},recruit_rate)
  //for h_exact=1
  out(0,5)=f_yi*f_di*gamma(0)*(di==di_rep);
  //for h_heap=1
  if(di_rep%5==0){
    out(1,5)=f_yi*f_di*gamma(1);
    if(recruit_rate>=0){//issued coupons->conditional reporting rules
      if(recruit_rate<ave_rate){ //recruit less than average
        if(di_rep>di){
          int k=(di_rep-5*floor(di/5))/5;
          out(1,5)=out(1,5)*pow(delta_1,k)/(1/(1-delta_1)-1);     
        }
        else //under assumption, di_rep>di if self recruitment rate<ave
          out(1,5)=0;
      }else{//recruit_rate>=ave_rate
        if(floor(di/5)>0){//di>=5
          double denom=1/(1-delta_2)-1 +//sum_n>0 delta_2^n 
            (1-pow(1-delta_2,floor(di/5)))/delta_2 ; // sum_{n>=0&<=K} (1-delta_2)^n
          int k=(di_rep-5*floor(di/5))/5;
          if(k>0)
            out(1,5)=out(1,5)*pow(delta_2,k)/denom;
          else
            out(1,5)=out(1,5)*pow(1-delta_2,-k)/denom;
        }
        else{//di<5
          double denom=1/(1-delta_2)-1 ; //sum_n>0 delta_2^n 
          int k=di_rep/5;
          out(1,5)=out(1,5)*pow(delta_2,k)/denom;
        }
      }
    }else{//recruit_rate=-1<0 wasn't issued coupons-> degenerating reporting rules, same logic as "recruit_rate>=ave_rate"
       if(floor(di/5)>0){//di>=5
          double denom=1/(1-delta_2)-1 +//sum_n>0 delta_2^n 
            (1-pow(1-delta_2,floor(di/5)))/delta_2 ; // sum_{n>=0&<=K} (1-delta_2)^n
          int k=(di_rep-5*floor(di/5))/5;
          if(k>0)
            out(1,5)=out(1,5)*pow(delta_2,k)/denom;
          else
            out(1,5)=out(1,5)*pow(1-delta_2,-k)/denom;
        }
        else{//di<5
          double denom=1/(1-delta_2)-1 ; //sum_n>0 delta_2^n 
          int k=di_rep/5;
          out(1,5)=out(1,5)*pow(delta_2,k)/denom;
        }
    }
  }
  else
    out(1,5)=0;

  //for h_error=1
  out(2,5)=f_yi*f_di*gamma(2);
  out(2,5)=out(2,5)*exp(lgamma(di_rep+phi)-lgamma(phi)-lgamma(di_rep+1)+
    di_rep*log(di/(di+phi))+phi*log(phi/(di+phi)));

  if(recruit_rate>=0){//issued coupons->conditional reporting rules
    if(recruit_rate<ave_rate){
      if(di_rep>di){
        //di-truncated Poisson (di)
        /*vector<double> num(di+1);
        iota (begin(num), end(num),0);
        for_each(num.begin(),num.end(),[&di](double& x){x=exp(-di+x*log(di)-logfactorial(x));});//Pr(x|di),0<=X<=di ~ Poisson(di)
        out(2,5)=out(2,5)*(di_rep>di)*exp(di_rep*log(di)-di-logfactorial(di_rep))/(1-accumulate(num.begin(), num.end(), 0.0));
        */
        //di-truncated NB(mean=di,phi)
        //NumericVector res=IbetaFromR(di/(di+phi),di+1,phi);
        //out(2,5)=out(2,5)/res[0]*R::beta(di+1,phi);
        
        out(2,5)=out(2,5)/boost::math::ibeta(di+1,phi,di/(di+phi));

        }
      else
        out(2,5)=0;
    }else {
      //0-truncated Poisson(di) 
      //out(2,5)=out(2,5)*(di_rep>0)*exp(di_rep*log(di)-di-logfactorial(di_rep))/(1-exp(-di));
      //0-truncated NB(mean=di,phi)
      out(2,5)=out(2,5)*(di_rep>0)/(1-pow(phi/(di+phi),phi));

    }
  }else{//recruit_rate=-1<0 wasn't issued coupons-> degenerating reporting rules, same logic as "recruit_rate>=ave_rate"
      out(2,5)=out(2,5)*(di_rep>0)/(1-pow(phi/(di+phi),phi));

  }
  return out; 
}

struct ExpDeg_i : public Worker
{   
  // input
  const vec draw_di;//all possible di's
  const int yi;
  const int di_rep;
  const double recruit_rate;
  const double ave_rate;
  const double pPAT;
  const arma::mat beta;
  const double mu;
  const double delta_1;
  const double delta_2;
  const double phi;
  // output: accumulated values  
  vec Emis_di;
  vec Emis_hi;
  vec mtx_beta_exact;
  vec mtx_beta_heap;
  arma::mat mtx_inverse_beta;
  double Emis_num_delta_1;
  double Emis_denom_delta_1;
  double Emis_S_delta_1;
  double Emis_S2_delta_1;
  double Emis_S_delta_2;
  double Emis_S2_delta_2;
  double Emis_S_phi;
  double Emis_S2_phi;
  // constructors
  // constructor1: The main constructor
  ExpDeg_i(const vec draw_di, const int yi,const int di_rep,
           const double recruit_rate, const double ave_rate,const double pPAT,const mat beta,
           const double mu,const double delta_1,const double delta_2,const double phi): 
    draw_di(draw_di), yi(yi), di_rep(di_rep),recruit_rate(recruit_rate),ave_rate(ave_rate),
    pPAT(pPAT), beta(beta),mu(mu),delta_1(delta_1),delta_2(delta_2),phi(phi),
    Emis_di(),Emis_hi(), mtx_beta_exact(),mtx_beta_heap(),mtx_inverse_beta(),
    Emis_num_delta_1(0.0),Emis_denom_delta_1(0.0),
    Emis_S_delta_1(0.0),Emis_S2_delta_1(0.0),Emis_S_delta_2(0.0),Emis_S2_delta_2(0.0),
    Emis_S_phi(0.0),Emis_S2_phi(0.0)
    {
    Emis_di.resize(2);Emis_di.fill(0.0);
    Emis_hi.resize(3);Emis_hi.fill(0.0);
    mtx_beta_exact.resize(2);mtx_beta_exact.zeros();
    mtx_beta_heap.resize(2);mtx_beta_heap.zeros();
    mtx_inverse_beta.resize(4,4);mtx_inverse_beta.zeros();
  }
  // constructor2: call for each split job
  ExpDeg_i(const ExpDeg_i& sum, Split) : 
    draw_di(sum.draw_di), yi(sum.yi), di_rep(sum.di_rep),
    recruit_rate(sum.recruit_rate),ave_rate(sum.ave_rate),
    pPAT(sum.pPAT), beta(sum.beta),mu(sum.mu),delta_1(sum.delta_1),delta_2(sum.delta_2),phi(sum.phi),
    Emis_di(),Emis_hi(),mtx_beta_exact(),mtx_beta_heap(), mtx_inverse_beta(),
    Emis_num_delta_1(0.0),Emis_denom_delta_1(0.0),
    Emis_S_delta_1(0.0),Emis_S2_delta_1(0.0),Emis_S_delta_2(0.0),Emis_S2_delta_2(0.0),
    Emis_S_phi(0.0),Emis_S2_phi(0.0)
    {
    Emis_di.resize(2);Emis_di.fill(0.0);
    Emis_hi.resize(3);Emis_hi.fill(0.0);
    mtx_beta_exact.resize(2);mtx_beta_exact.zeros();
    mtx_beta_heap.resize(2);mtx_beta_heap.zeros();
    mtx_inverse_beta.resize(4,4);mtx_inverse_beta.zeros();
  }
  // accumulate just the element of the range I've been asked to
  void operator() (size_t begin,size_t end) {
    for (size_t i = begin; i < end; i++){
      int di=draw_di(static_cast<R_xlen_t>(i));
      mat prob=cal_complete(di,mu,beta,yi,di_rep,recruit_rate,ave_rate,pPAT,delta_1,delta_2,phi);
      Emis_di[0] += accu(prob.col(5));//sum_j Pr(y_obs,y_mis={di,h_{ij}=1})=Pr(y_obs,y_mis={di})
      Emis_di[1] += di*accu(prob.col(5));
      Emis_hi[0] += prob(0,5);//Pr(y_obs,y_mis={di,h_{i,exact}=1})
      Emis_hi[1] += prob(1,5);//Pr(y_obs,y_mis={di,h_{i,heap}=1})
      Emis_hi[2] += prob(2,5);//Pr(y_obs,y_mis={di,h_{i,error}=1})
      
      colvec Ai={1,di};
      mtx_beta_exact += prob(0,5)*(1-prob(0,4))*Ai+//(1-gamma_exact)Ai*Pr(y_obs,y_mis={di,h_{i,exact}=1})
        (prob(1,5)+prob(2,5))*(-prob(0,4))*Ai;//(-gamma_exact)Ai*Pr(y_obs,y_mis={di,h_{i,exact}=0})
      mtx_beta_heap += prob(1,5)*(1-prob(1,4))*Ai+(prob(0,5)+prob(2,5))*(-prob(1,4))*Ai;
      mat add(4,4);
      add.submat(0,0,1,1)= -prob(0,4)*(1-prob(0,4))*Ai*Ai.t();
      add.submat(0,2,1,3)= prob(0,4)*prob(1,4)*Ai*Ai.t();
      add.submat(2,2,3,3)= -prob(1,4)*(1-prob(1,4))*Ai*Ai.t();
      add.submat(2,0,3,1)= prob(0,4)*prob(1,4)*Ai*Ai.t();
      mtx_inverse_beta += add*accu(prob.col(5));
      
      if(di_rep%5==0 & di_rep>di & recruit_rate<ave_rate){//conditions for delta_1
        int k=(di_rep-5*floor(di/5))/5;
        Emis_num_delta_1+=prob(1,5);
        Emis_denom_delta_1+=k*prob(1,5);//*Pr(y_obs,y_mis={di,h_{i,heap}=1})
        Emis_S_delta_1+=(-1/(1-delta_1)+(k-1)/delta_1)*prob(1,5);//*Pr(y_obs,y_mis={di,h_{i,heap}=1})
        Emis_S2_delta_1+=(-1/pow(1-delta_1,2)-(k-1)/pow(delta_1,2))*prob(1,5);
      }
      if(di_rep%5==0 & recruit_rate>=ave_rate){//conditions for delta_2
        int k=(di_rep-5*floor(di/5))/5;
        if(di<5){
          Emis_S_delta_2+=(-1/(1-delta_2)+(k-1)/delta_2)*prob(1,5);//*Pr(y_obs,y_mis={di,h_{i,heap}=0})
          Emis_S2_delta_2+=(-1/pow(1-delta_2,2)-(k-1)/pow(delta_2,2))*prob(1,5);
        }
        if(di>=5){
          int K=floor(di/5)-1;
          double temp1=2*delta_2-1+pow(1-delta_2,K+2)*(K*delta_2+1);
          double temp2=1+(delta_2-1)*(delta_2+pow(1-delta_2,K+1));
          double part_S=-temp1/((1-delta_2)*delta_2*temp2);
          double part_S2=-2/((1-delta_2)*delta_2*temp2)+
            (pow(K,2)*delta_2+3*K*delta_2+2)*pow(1-delta_2,K)/(delta_2*temp2)+
            temp1*(1-2*delta_2)/pow(delta_2,2)/pow(1-delta_2,2)/temp2+
            temp1*(pow(1-delta_2,K)*(K+2-2*delta_2-delta_2*K)+2*delta_2-1)/(1-delta_2)/delta_2/pow(temp2,2);
          if(k>0){
            Emis_S_delta_2+=(k/delta_2+part_S)*prob(1,5);
            Emis_S2_delta_2+=(-k/pow(delta_2,2)+part_S2)*prob(1,5);  		
          }
          else{
            k=-k;
            Emis_S_delta_2+=(-k/(1-delta_2)+part_S)*prob(1,5);
            Emis_S2_delta_2+=(-k/pow(1-delta_2,2)+part_S2)*prob(1,5);
          }
        }
      }
      
      if(recruit_rate>=ave_rate&di_rep>0){//conditions for 0-truncated NB 
        vec num = linspace<vec>(0, di_rep-1,di_rep);
        num.transform([&](double val){return (1/(val+phi)); } );
        double temp1=1/(1-pow(phi/(di+phi),phi))-1;
        double temp2=di/(di+phi)+log(phi/(di+phi));
        Emis_S_phi+=(accu(num)+(di-di_rep)/(di+phi)+log(phi/(di+phi))+
          temp1*temp2)*prob(2,5);//*Pr(y_obs,y_mis={di,h_{i,error}=1})
        num = linspace<vec>(0, di_rep-1,di_rep);
        num.transform([&](double val){return (-1/pow(val+phi,2)); } );
        Emis_S2_phi+=(accu(num)+(di_rep-di)/pow(di+phi,2)+di/phi/(di+phi)+
          temp1*pow(di,2)/phi/pow(di+phi,2)+
          pow(phi/(di+phi),phi)*pow(temp2,2)/pow(1-pow(phi/(di+phi),phi),2))*prob(2,5);
       }
      if(recruit_rate<ave_rate&di_rep>di){//conditions for di-truncated NB 
        vec num = linspace<vec>(0, di_rep-1,di_rep);
        num.transform([&](double val){return (1/(val+phi)); } );
        vec num2 = linspace<vec>(0, di,di+1);
        num2.transform([&](double val){return (1/(val+phi)); } );
        double temp1=exp((phi-1)*log(phi)+(di+1)*log(di)-(di+phi+1)*log(di+phi));
        double temp3=boost::math::beta(di+1,phi,di/(di+phi));
        double temp2=IntegralFunc1(di/(di+phi),di,phi);
        double temp4=IntegralFunc2(di/(di+phi),di,phi);
        Emis_S_phi+=(accu(num)-accu(num2)+(di-di_rep)/(di+phi)+log(phi/(di+phi))-
          (-temp1+temp2)/temp3)*prob(2,5);
        
        num = linspace<vec>(0, di_rep-1,di_rep);
        num.transform([&](double val){return (-1/pow(val+phi,2)); } );
        num2 = linspace<vec>(0, di,di+1);
        num2.transform([&](double val){return (-1/pow(val+phi,2)); } );
        Emis_S2_phi+=(accu(num)-accu(num2)+(di_rep-di)/pow(di+phi,2)+di/phi/(di+phi)+
          pow((-temp1+temp2)/temp3,2)-(temp4-temp1*(2*log(phi/(di+phi))-1/phi-1/(di+phi)))/temp3
                        )*prob(2,5);
       
       }
     
    }
  }
  
  
  // join my value with that of another ExpDeg
  void join(const ExpDeg_i& rhs) {
    Emis_di[0]+=rhs.Emis_di[0];
    Emis_di[1]+=rhs.Emis_di[1];
    Emis_hi[0]+=rhs.Emis_hi[0];
    Emis_hi[1]+=rhs.Emis_hi[1];
    Emis_hi[2]+=rhs.Emis_hi[2];
    
    mtx_beta_exact += rhs.mtx_beta_exact;
    mtx_beta_heap += rhs.mtx_beta_heap;
    mtx_inverse_beta += rhs.mtx_inverse_beta;
    
    Emis_num_delta_1+= rhs.Emis_num_delta_1;
    Emis_denom_delta_1+= rhs.Emis_denom_delta_1;
    Emis_S_delta_1+= rhs.Emis_S_delta_1;
    Emis_S2_delta_1+= rhs.Emis_S2_delta_1;
    Emis_S_delta_2+= rhs.Emis_S_delta_2;
    Emis_S2_delta_2+= rhs.Emis_S2_delta_2;
    Emis_S_phi+= rhs.Emis_S_phi;
    Emis_S2_phi+= rhs.Emis_S2_phi;
  }
};

struct ExpDegList : public Worker
{   
  // input
  const vec YMale;
  const vec RepDeg;
  const arma::mat X;
  const vec alpha;
  const mat beta;
  const double pPAT;
  const vec draw_di;
  const vec recruit_rate_list;
  const double ave_rate;
  const double delta_1;
  const double delta_2;
  const double phi;
  // output: accumulated values numerator and denominator 
  double PrObs;
  vec Emis_di_list;
  mat Emis_hi_list;
  vec mtx_beta_exact;
  vec mtx_beta_heap;
  mat mtx_inverse_beta;
  double Emis_num_delta_1;
  double Emis_denom_delta_1;
  double Emis_S_delta_1;
  double Emis_S2_delta_1;
  double Emis_S_delta_2;
  double Emis_S2_delta_2;
  double Emis_S_phi;
  double Emis_S2_phi;
  // constructor: the RMatrix class can be automatically converted to from the Rcpp matrix type
  ExpDegList(const vec YMale, const vec RepDeg, const mat X,
             const vec alpha, const mat beta, const double pPAT,
             const vec draw_di, const vec recruit_rate_list,const double ave_rate,
             const double delta_1,const double delta_2,const double phi): 
    YMale(YMale), RepDeg(RepDeg), X(X), alpha(alpha), beta(beta), pPAT(pPAT), draw_di(draw_di),
    recruit_rate_list(recruit_rate_list),ave_rate(ave_rate),
    delta_1(delta_1), delta_2( delta_2),phi(phi),
    PrObs(),Emis_di_list(),Emis_hi_list(),mtx_beta_exact(),mtx_beta_heap(),
    mtx_inverse_beta(), Emis_num_delta_1(0.0),Emis_denom_delta_1(0.0),
    Emis_S_delta_1(0.0),Emis_S2_delta_1(0.0),Emis_S_delta_2(0.0),Emis_S2_delta_2(0.0),
    Emis_S_phi(0.0),Emis_S2_phi(0.0)
  {
    Emis_di_list.resize(X.n_rows);Emis_di_list.fill(0.0);
    Emis_hi_list.resize(X.n_rows,3);Emis_hi_list.fill(0.0);
    mtx_beta_exact.resize(2);mtx_beta_exact.zeros();
    mtx_beta_heap.resize(2);mtx_beta_heap.zeros();
    mtx_inverse_beta.resize(4,4);mtx_inverse_beta.zeros();
  };
  ExpDegList(const ExpDegList& sum, Split) :
    YMale(sum.YMale), RepDeg(sum.RepDeg), X(sum.X), 
    alpha(sum.alpha), beta(sum.beta), pPAT(sum.pPAT), draw_di(sum.draw_di),
    recruit_rate_list(sum.recruit_rate_list),ave_rate(sum.ave_rate),
    delta_1(sum.delta_1),delta_2(sum.delta_2),phi(sum.phi),
    PrObs(),Emis_di_list(),Emis_hi_list(), mtx_beta_exact(),mtx_beta_heap(),
    mtx_inverse_beta(), Emis_num_delta_1(0.0),Emis_denom_delta_1(0.0),
    Emis_S_delta_1(0.0),Emis_S2_delta_1(0.0),Emis_S_delta_2(0.0),Emis_S2_delta_2(0.0),
    Emis_S_phi(0.0),Emis_S2_phi(0.0)
  {
    Emis_di_list.resize(sum.X.n_rows);Emis_di_list.fill(0.0);
    Emis_hi_list.resize(sum.X.n_rows,3);Emis_hi_list.fill(0.0);
    mtx_beta_exact.resize(2);mtx_beta_exact.zeros();
    mtx_beta_heap.resize(2);mtx_beta_heap.zeros();
    mtx_inverse_beta.resize(4,4);mtx_inverse_beta.zeros();
  };
  
  void operator() (size_t begin,size_t end) {
    
    for (size_t i = begin; i < end; i++){
      int yi=YMale(static_cast<R_xlen_t>(i));
      int di_rep=RepDeg(static_cast<R_xlen_t>(i));
      double recruit_rate=recruit_rate_list(static_cast<R_xlen_t>(i));
      //double ave_rate=sum(recruit_rate_list)/(double)recruit_rate_list.size();
      rowvec xi=X.row(static_cast<R_xlen_t>(i));
      double mu = exp(inner_product(alpha.begin(),alpha.end(),xi.begin(),0.0));

      ExpDeg_i *di=new ExpDeg_i(draw_di,yi,di_rep,recruit_rate,ave_rate,pPAT,beta,mu,delta_1,delta_2,phi); 
      parallelReduce(0, draw_di.size(), *di);
       if(di->Emis_di[0]!=0){
        PrObs+=di->Emis_di[0];
        Emis_di_list[i]=di->Emis_di[1]/di->Emis_di[0];
        Emis_hi_list.row(i)=trans(di->Emis_hi/di->Emis_di[0]);
        mtx_beta_exact += di->mtx_beta_exact/di->Emis_di[0];
        mtx_beta_heap += di->mtx_beta_heap/di->Emis_di[0];
        mtx_inverse_beta += di->mtx_inverse_beta/di->Emis_di[0];
        Emis_num_delta_1 += di->Emis_num_delta_1/di->Emis_di[0]; 
        Emis_denom_delta_1 += di->Emis_denom_delta_1/di->Emis_di[0]; 
        Emis_S_delta_1 += di->Emis_S_delta_1/di->Emis_di[0]; 
        Emis_S2_delta_1 += di->Emis_S2_delta_1/di->Emis_di[0];
        Emis_S_delta_2 += di->Emis_S_delta_2/di->Emis_di[0]; 
        Emis_S2_delta_2 += di->Emis_S2_delta_2/di->Emis_di[0]; 
        Emis_S_phi += di->Emis_S_phi/di->Emis_di[0]; 
        Emis_S2_phi += di->Emis_S2_phi/di->Emis_di[0]; 
      }else{
        Emis_di_list[i]=0;
      }
      delete di;
      di=NULL;
    }
  }
  void join(const  ExpDegList& rhs) {
    PrObs+=rhs.PrObs;
    Emis_di_list+=rhs.Emis_di_list;
    Emis_hi_list+=rhs.Emis_hi_list;
    mtx_beta_exact+=rhs.mtx_beta_exact;
    mtx_beta_heap+=rhs.mtx_beta_heap;
    mtx_inverse_beta+=rhs.mtx_inverse_beta;
    Emis_num_delta_1+= rhs.Emis_num_delta_1;
    Emis_denom_delta_1+= rhs.Emis_denom_delta_1;
    Emis_S_delta_1+= rhs.Emis_S_delta_1;
    Emis_S2_delta_1+= rhs.Emis_S2_delta_1;
    Emis_S_delta_2+= rhs.Emis_S_delta_2;
    Emis_S2_delta_2+= rhs.Emis_S2_delta_2;
    Emis_S_phi+= rhs.Emis_S_phi;
    Emis_S2_phi+= rhs.Emis_S2_phi;
  }
};

struct update_ALPHA: public Worker
{   
  // input
  const vec old_alpha;
  const vec YMale;
  const mat X;
  const vec Emis_di_list;
  // output: accumulated values numerator&denominator 
  mat mtx_inverse_alpha;
  vec mtx_alpha;
  
  // constructors
  // constructor1: The main constructor
  update_ALPHA(const vec old_alpha,
               const vec YMale, const mat X, const vec Emis_di_list): 
    old_alpha(old_alpha), 
    YMale(YMale), X(X), Emis_di_list(Emis_di_list),
    mtx_inverse_alpha(),mtx_alpha()
  {
    mtx_inverse_alpha.resize(X.n_cols,X.n_cols);mtx_inverse_alpha.zeros();
    mtx_alpha.resize(X.n_cols);mtx_alpha.zeros();
  }
  // constructor2: call for each split job
  update_ALPHA(const update_ALPHA& sum, Split) : 
    old_alpha(sum.old_alpha), 
    YMale(sum.YMale), X(sum.X), Emis_di_list(sum.Emis_di_list),
    mtx_inverse_alpha(),mtx_alpha()
  {
    mtx_inverse_alpha.resize(sum.X.n_cols,sum.X.n_cols);mtx_inverse_alpha.zeros();
    mtx_alpha.resize(sum.X.n_cols);mtx_alpha.zeros();
  }
  
  // accumulate just the element of the range I've been asked to
  void operator() (size_t begin,size_t end) {
    
    for (size_t i = begin; i < end; i++){
      rowvec xi=X.row(static_cast<R_xlen_t>(i));
      double mu = exp(inner_product(old_alpha.begin(),old_alpha.end(),xi.begin(),0.0));
      mtx_inverse_alpha += mu*(exp(-mu)*(3-mu-2*exp(-mu))-1)/pow(1-exp(-mu),2)*xi.t()*xi;
      mtx_alpha += (Emis_di_list[static_cast<R_xlen_t>(i)]-2*mu+mu/(1-exp(-mu)))*xi.t();
    }
  }
  void join(const update_ALPHA& rhs) {
    mtx_inverse_alpha += rhs.mtx_inverse_alpha;
    mtx_alpha += rhs.mtx_alpha;
  }
};

double max_diff_param(vec param_old_alpha,vec param_new_alpha,
                      mat param_old_beta,mat param_new_beta,
                      double delta_1_old, double delta_1,
                      double delta_2_old, double delta_2, double phi_old,double phi
                        ){  
  vec alpha_dist=param_old_alpha-param_new_alpha;
  alpha_dist.transform( [](double val) { return (val<0 ? -val : val); } );
  
  mat beta_dist=param_old_beta-param_new_beta;
  beta_dist.transform( [](double val) { return (val<0 ? -val : val); } );
  
  return max(
           max(
            max(
              max(alpha_dist.max(),beta_dist.max()),
              abs(delta_1_old-delta_1) ),
            abs(delta_2_old-delta_2) ),
           abs(phi_old-phi))
           ;
}

//[[Rcpp::export]]
List runEM(const NumericVector draw_di0,
           const int itermax,
           const double restrict,
           const double threshold,
           const NumericVector alpha0,
           const NumericMatrix beta0,
           const double pPAT,
           const NumericMatrix X0,
           const NumericVector YMale0,
           const NumericVector RepDeg0,
           const NumericVector recruit_rate_list0,
           const double ave_rate,
           double delta_1,double delta_2,double phi){

  //type tranfer
  vec draw_di = as<vec>(draw_di0);
  vec alpha = as<vec>(alpha0);
  mat beta = as<mat>(beta0);
  mat X = as<mat>(X0);
  vec YMale = as<vec>(YMale0);
  vec RepDeg = as<vec>(RepDeg0);
  vec recruit_rate_list = as<vec>(recruit_rate_list0);
  
  int iter=0; 
  list<double> di_diff_iter;
  list<double> PrObs_iter;
  list<stdvec> alpha_iter;
  list<stdvecvec> beta_iter;
  list<stdvec> Emis_di_list_iter;
  list<stdvecvec>Emis_hi_list_iter;
  list<double> delta_1_iter;list<double> delta_2_iter;
  list<double> phi_iter;
  alpha_iter.push_back(conv_to< stdvec >::from(alpha)); 
  beta_iter.push_back(mat_to_std_vec(beta));
  delta_1_iter.push_back(delta_1);delta_2_iter.push_back(delta_2);
  phi_iter.push_back(phi);
  
  double diff=1;
  vec di_old=RepDeg;
  vec temp=RepDeg;
  double di_diff=0;
  vec alpha_old=alpha;
  mat beta_old=beta;
  double delta_1_old=delta_1; double delta_2_old=delta_2;
  double phi_old=phi;
  
  vector<double> not_converge(6,0);
  while(iter<itermax&&diff>threshold){
    ExpDegList *expdeglist= new ExpDegList(YMale,RepDeg,X,alpha,beta,pPAT,draw_di,recruit_rate_list,ave_rate,delta_1,delta_2,phi);
    parallelReduce(0,X.n_rows, *expdeglist);
    PrObs_iter.push_back(expdeglist->PrObs); 
    delta_1 =1- expdeglist->Emis_num_delta_1/expdeglist->Emis_denom_delta_1;
    if(delta_1<0|delta_1>1)delta_1=0.5;
    /*double Emis_S_delta_1 = expdeglist->Emis_S_delta_1;
    double Emis_S2_delta_1 = expdeglist->Emis_S2_delta_1;
    delta_1 = delta_1_old-Emis_S_delta_1/Emis_S2_delta_1;*/
    delta_1_iter.push_back(delta_1);
    double Emis_S_delta_2 = expdeglist->Emis_S_delta_2;
    double Emis_S2_delta_2 = expdeglist->Emis_S2_delta_2;
    delta_2 = delta_2_old-1/Emis_S2_delta_2*Emis_S_delta_2;
    if(delta_2<0|delta_2>1)delta_2=0.5;
    delta_2_iter.push_back(delta_2);
    
    double Emis_S_phi = expdeglist->Emis_S_phi;
    double Emis_S2_phi = expdeglist->Emis_S2_phi;
    phi = phi_old-1/Emis_S2_phi*Emis_S_phi;
    if(phi>1e10)phi=1e10;
    if(phi<0)phi=1;
    phi_iter.push_back(phi);
    
    vec mtx_beta = join_cols(expdeglist->mtx_beta_exact,expdeglist->mtx_beta_heap);
    vec mtx_beta_prod = -inv(expdeglist->mtx_inverse_beta)*mtx_beta;
    beta.row(0)=beta.row(0) + mtx_beta_prod.subvec(0,1).t();
    beta.row(1)=beta.row(1) + mtx_beta_prod.subvec(2,3).t();
    beta.transform( [&restrict ](double val) { return (val<-restrict ? -restrict : val); } );
    beta.transform( [&restrict ](double val) { return (val>restrict ? restrict : val); } );
    
    update_ALPHA *update_alpha = new update_ALPHA(alpha,YMale,X,expdeglist->Emis_di_list);
    parallelReduce(0, X.n_rows, *update_alpha);
    alpha=alpha-inv(update_alpha->mtx_inverse_alpha)*update_alpha->mtx_alpha;
 
    delete update_alpha;
    update_alpha=NULL;  
    
    alpha_iter.push_back(conv_to< stdvec >::from(alpha)); 
    beta_iter.push_back(mat_to_std_vec(beta));
    Emis_di_list_iter.push_back(conv_to< stdvec >::from(expdeglist->Emis_di_list));
    Emis_hi_list_iter.push_back(mat_to_std_vec(expdeglist->Emis_hi_list));
    
    temp=expdeglist->Emis_di_list;
    temp=di_old-temp;
    di_diff=dot(temp.t(),temp);
    di_diff=pow(di_diff,.5)/pow(dot(di_old.t(),di_old),.5);
    di_diff_iter.push_back(di_diff);
    di_old=expdeglist->Emis_di_list;
    
    
    
    
    iter=iter+1;
    diff=max_diff_param(alpha_old,alpha,beta_old,beta,
                        delta_1_old,delta_1,delta_2_old,delta_2,phi_old,phi
                          );
    if(iter==itermax&&diff>threshold){//not converge until last iteration
      
      vec alpha_dist=alpha_old-alpha;
      alpha_dist.transform( [](double val) { return (val<0 ? -val : val); } );
      if(alpha_dist.max()>threshold)
        not_converge[0]=alpha_dist.max();
      mat beta_dist=beta_old-beta;
      beta_dist.transform( [](double val) { return (val<0 ? -val : val); } );
      if(beta_dist.max()>threshold)
        not_converge[1]=beta_dist.max();
      if(abs(delta_1_old-delta_1)>threshold)
        not_converge[2]=abs(delta_1_old-delta_1);
      if(abs(delta_2_old-delta_2)>threshold)
        not_converge[3]=abs(delta_2_old-delta_2);
      if(abs(phi_old-phi)>threshold)
        not_converge[4]=abs(phi_old-phi);
    }
    
    
    
    alpha_old=alpha;
    beta_old=beta;
    delta_1_old=delta_1; 
    delta_2_old=delta_2;
    phi_old=phi;
    
  }
  
  list<double>::iterator it1 = PrObs_iter.begin();
  list<double>::iterator it2 = PrObs_iter.begin();
  // Move the iterator 
  advance(it1, PrObs_iter.size()-1);
  advance(it2, PrObs_iter.size()-2);
  // Print the element at the it
  not_converge[5]=*it1-*it2;
  
  List out;
  out=List::create(Named("degree")=di_old,
                   _["alpha"]=alpha,
                   _["beta"]=beta,
                   _["delta_1"]=delta_1,
                   _["delta_2"]=delta_2,
                   _["phi"]=phi);
  return out;
}

