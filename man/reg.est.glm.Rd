% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_regression.R
\name{reg.est.glm}
\alias{reg.est.glm}
\title{Logistic regression modeling in a RDS study}
\usage{
reg.est.glm(data, formula, include.degree = F)
}
\arguments{
\item{data}{Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest}

\item{formula}{An object of class 'formula'; description of the model with dependent and independent variables. Note that for linear regression the dependent variable must be continuous. Use reg.est.lm to perform regression with continuous dependent variable.}

\item{include.degree}{Logical; when TRUE the function performs weighted logistic regression with DEGREE used as weights.}
}
\value{
A data frame consisting of the point estimates and associated standard errors
}
\description{
Logistic regression modeling in a RDS study
}
\examples{
reg_glm_out = reg.est.glm(data, formula, include.degree=F)
}
