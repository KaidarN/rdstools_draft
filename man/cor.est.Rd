\name{correlation_pearson}
\alias{correlation_pearson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation analysis in a RDS study
}
\description{
This function performs Pearson correlation on the bootstraped data.
}
\usage{
out = cor.est(x, y, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector; a continuous variable x}
  \item{y}{Numeric vector; a continuous variable y}
  \item{data}{Data frame; Contains data from bootstrap_RDS function with boot_n and variables of interest}
}
}
\value{
A numeric scalar value indicating correlation coefficient between x and y
}
\author{
THE AUTHORS WILL BE PLACED HERE
}

\examples{
EXAMPLES WILL BE PLACED HERE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
