\name{linear_regression}
\alias{linear_regression}

\title{
Linear regression modeling in a RDS study
}
\description{
This function performs weighted and unweighted multiple linear regression modeling on the bootstraped data.
}
\usage{
out = reg.est.lm(data, formula, include.degree = F)
}

\arguments{
  \item{data}{Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest
}
  \item{formula}{An object of class 'formula'; description of the model with dependent and independent variables. Note that for linear regression the dependent variable must be continuous. Use reg.est.glm to perform regression with binary categorical dependent variable.}
  \item{include.degree}{Logical; when TRUE the function performs weighted linear regression with DEGREE used as weights}
}

\value{
A data frame consisting of the point estimates and associated standard errors
}

\author{
THE AUTHORS WILL BE PLACED HERE
}

\examples{
EXAMPLES WILL BE PLACED HERE
}
