\name{plot_beta_curve}
\alias{plot_beta_curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate power curves for beta-distributed case/control measurements
}
\description{plot_beta_curve() plots power estimates on the y axis and effect size on the x axis for each of the case:control ratios tested.
}
\usage{
plot_beta_curve(data, a, MHT, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe produced with beta_curve()
}
\item{a}{the p-value significance level used to generate the data
}
\item{MHT}{the multiple hypotheses testing correction method used to generate the data
}
\item{q}{the adjusted p-value significance level used to generate the data
}
}
\examples{
#Generate power curves for beta-distributed case/control measurements
my_beta_data <- beta_curve(mu0=0.2, sd0=0.15, mu1_min=0.25,
mu1_max=0.27, mu1_step=0.01, sd1=0.15, n0s=c(150, 150, 150),
n1s=c(150,300, 450),N=100,a=0.05,MHT="bonferroni",q=0.1)
plot_beta_curve(data=my_beta_data, a=0.05, MHT="bonferroni",q=0.1)
}
