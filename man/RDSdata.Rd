\name{RDS_data_handling_weight}
\alias{RDS_data_handling_weight}

\title{
Data processing in a RDS study
}
\description{
This function is used to process the data from an RDS study. The function extracts unique_id, redeemed coupon number and issued coupons number from the original data. Processing data with RDSdata function is neccessary for appropopriate work of bootstrap function.
}
\usage{
out = RDSdata(data, unique_id, redeemed_coupon, issued_coupon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame; Should contain ID numbers for the nodes in the social network, corresponding redeemed coupon number, and issued coupon number}
  \item{unique_id}{Character vector; The column name of the column in the data that represents the ID numbers for the nodes in the social network}
  \item{redeemed_coupon}{Character vector; The column name of the column in the data that represents coupon numbers of coupons redeemed by respondents when participating in the survey}
  \item{issued_coupon}{Character vector; The column name of the column in the data that represents the coupon numbers of coupons issued to respondents}
  \item{degree}{Numeric vector;  The column name of the column in the data that represents the degree of the respondents}
  \item{degree_miss}{Character; This parameter is used to set the method for imputing missing values in the 'degree' variable. There are three methods to choose from: mean, median, and hotdeck. If this parameter is not set, the default method is hotdeck.}
  \item{degree_zero}{Character; This parameter is used to set the method for imputing zero values in the 'degree' variable. Three methods are available for selection: mean, median, and hotdeck. If this parameter is not set, the default imputation method is hotdeck}
}
\details{
RDS data function vingette provide a more detailed explanation of how to use the data preprocessing.
}
\value{
A data frame consisting of the following elements:
  \item{ID}{Character vector; Renamed unique_id variable}
  \item{R_CP}{Character vector; Renamed redeemed coupon variable}
  \item{T_CP1 - T_CPn}{Character vector(s); Renamed issued coupon variable}
  \item{DEGREE}{Numeric vector; Renamed degree variable}
  \item{WAVE}{Numeric vector; Indicates which round the node was introduced into the survey, and the value of Seed is 0}
  \item{S_ID}{chacarter vector; Indicates the ID of the seed corresponding to the node. The value of the seed is itself}
  \item{R_ID}{Character vector; Indicates the ID of the node who recruits the node joinning the survey. The value of seed is NA}
  \item{SEED}{Numeric vector: Values are only 0 and 1, they are used to indicate whether the node is seed or not. If it is seed, the value is 1, if not, it is 0.}
  \item{CT_T_CT}{Numeric data, indicates how many coupons have been issued to this node to invite others to join the survey}
  \item{CT_T_CP_USED}{Numeric data, indicates how many coupons issued to this node have been used to invite others to join the survey.}
  \item{DEGREE_MODIFY}{Numeric data, indicates the `degree` data after imputation of zero values and missing values.}
}

\author{
THE AUTHORS WILL BE PLACED HERE
}


\examples{
EXAMPLES WILL BE PLACED HERE
}

