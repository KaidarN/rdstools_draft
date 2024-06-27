\name{Bootstrap tree and chain functions}
\alias{Bootstrap tree and chain functions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance estimation with bootstrap chain and tree methods
}
\description{
This function implements bootstrap chain, unidirectional bootstrap tree and bidirectional bootstrap tree methods for variance estimation in an RDS study
}
\usage{
out = bootstrap_RDS(RESPONDENT_ID, SEED_ID, SEED, RECRUITER_ID,
                         data.cov, type, n.times, return_data=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RESPONDENT_ID}{character or numeric vector; A variable indicating respondent ID}
  \item{SEED_ID}{character or numeric vector; A variable indicating seed ID}
  \item{RECRUITER_ID}{character or numeric vector; A variable indicating recruiter ID}
  \item{data.cov}{vector, matrix or data frame; Variable(s) of interest that will be used in the analysis}
  \item{type}{character; One of the six types of bootstrap methods, please see the vignette for a more detailed information}
}
  \item{n.times}{scalar; How many bootstrap resamples should be taken}
  \item{return_data}{logical; By default return_data=F returns respondent ID and bootstrap iterations. If return_data TRUE, then the output will include variable(s) specified in data.cov}
}
\details{
Bootstrap vingettes provide more detailed explanation for each bootstrap method implemented in bootstrap_RDS
}
\value{
When return_data=F returns a data frame consistiting of the following elements:
  \item{RESPONDENT_ID}{character vector; A variable indicating respondent ID}
  \item{boot_n}{character vector; A indicator variable for each bootstrap sample}
}

\examples{
Some examples here
}

