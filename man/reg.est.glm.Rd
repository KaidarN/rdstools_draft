\name{logistic_regression}
\alias{logistic_regression}

\title{
Logistic regression modeling in a RDS study
}
\description{
This function performs weighted and unweighted multiple logistic regression modeling on the bootstraped data.
}
\usage{
out = reg.est.glm(data, formula, include.degree = F)
}

\arguments{
  \item{data}{Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest}
  \item{formula}{An object of class 'formula'; description of the model with dependent and independent variables. Note that for linear regression the dependent variable must be continuous. Use reg.est.lm to perform regression with continuous dependent variable.}
  \item{include.degree}{Logical; when TRUE the function performs weighted logistic regression with DEGREE used as weights.}
}

\value{
A data frame consisting of the point estimates and associated standard errors
}

\author{
THE AUTHORS WILL BE PLACED HERE
}

\examples{
THE EXAMPLES WILL BE PLACED HERE
}

