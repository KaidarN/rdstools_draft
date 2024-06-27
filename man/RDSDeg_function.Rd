\name{RDSDeg_Function}
\alias{RDSDeg_Function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating individual degree size in a RDS study
}
\description{
 DESCRIPTION WILL BE HERE
}
\usage{
out <- Est_Deg(sample,pPAT=0.01,
           covariates=c("Male","Black","Married","LowEduc","Age","Income","DETROIT",
                         "M_DQ3","M_DQ7","M_DQ29_SUPPORT","M_GQ11","M_HQ9","Duration"),
             maxiter=1000,beta_restrict=10000,threshold=0.005)
}

\arguments{
  \item{sample}{data.frame; An RDS sample, including information of: ID (Respondent unique ID), recruiterID (Recruiter's ID, not necessary for this function), seedID (Seed ID, not necessary for this function), wave (Wave number, not necessary for this function), RepDeg (Respondent's reported degree), YPat (Respondent's reported peers named Pat within target RDS population ), RecRate (Respondent's recruitment rate, i.e., the ratio of the number of recruits to the number of issued coupons ), Covariates (Subject-specific characteristics. If there exists continuous variables with large scales, we recommend a log- transform).
}
  \item{pPat}{scalar; External knowledge of the  size proportion of subpopulation named "Pat" within the target RDS population.}
  \item{covariates}{vector; A vector of subject-specific characteristics that are assumed to be associated with the true degree distribution via a Poisson regression model.}
  \item{maxiter}{scalar; Maximum number of iterations used in EM algorithm. Default is 1000.}
  \item{beta_restrict}{scalar; Bound of coefficients from the parametric regression model of the underlying reporting behavior. Default is 10000.}
  \item{threshold}{scalar; Convergence threshold. Default is 0.005.}


}
\details{
  @useDynLib RDStools, .registration = TRUE
  @importFrom Rcpp sourceCpp
}
\value{
A list consisting of the following elements:
  \item{degree}{vector; estimated individual true degree size.}
  \item{alpha}{vector; Estimated coefficients obtained from the parametric regression model of the true degree adjusted for subject-specific characteristics.}
  \item{beta}{matrix; Estimated coefficients obtained from the parametric regression model of the underlying reporting behavior adjusted for the individual true degree. }
  \item{delta_1, delta_2, phi}{scalars; Estimated parameters obtained from models of reporting rules.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Yibo Wang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
out <- Est_Deg(sample,pPAT=0.01,
               covariates=c("Male","Black","Married","LowEduc","Age","Income","DETROIT",
               "M_DQ3","M_DQ7","M_DQ29_SUPPORT","M_GQ11","M_HQ9","Duration"),
               maxiter=1000,beta_restrict=10000,threshold=0.005)
}

