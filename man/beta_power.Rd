\name{beta_power}
\alias{beta_power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate corrected power to detect a certain effect size between beta-distributed cases and controls at a certain significance level}
\description{beta_power() calculates power as the fraction of p-values or adjusted p-values of the beta regression's case/control coefficient in N simulations of case/control beta numbers below a predefined threshold. Allows for case and control groups of unequal size. Outputs a list of the power estimate based on the p-value, a vector of all p-values, the power estimate based on the adjusted p-value and a vector of all adjusted p-values.
}
\usage{
beta_power(mu0, sd0, mu1, sd1, n0, n1, N,a,MHT,q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{a reasonable estimate of the control group's mean
}
\item{sd0}{a reasonable estimate of the control group's standard deviation
}
\item{mu1}{a reasonable estimate of the case group's mean
}
\item{sd1}{a reasonable estimate of the case group's standard deviation
}
\item{n0}{the sample size of the control group
}
\item{n1}{the sample size of the case group
}
\item{N}{the number of simulations to be performed (e.g. 1,000)
}
\item{a}{the significance level at which the p-value will be tested
}
\item{MHT}{the correction method used for multiple hypotheses testing (e.g. "BH","bonferroni"). For all methods available: ?p.adjust.
}
\item{q}{the significance level at which the adjusted p-value will be tested
}
}
\examples{
#Calculate power for a certain effect size
my_power <- beta_power(mu0=0.2, sd0=0.15, mu1=0.25, sd1=0.15,
n0=150, n1=300, N=100,a=0.05,MHT="bonferroni",q=0.1)
p_value_power <- my_power[[1]]
p_values <- my_power[[2]]
adjusted_p_value_power <- my_power[[3]]
adjusted_p_values <- my_power[[4]]
}
