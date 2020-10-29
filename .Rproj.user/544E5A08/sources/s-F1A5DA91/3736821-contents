\name{mu1_curve_data}
\alias{mu1_curve_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate power to detect a range of effect sizes between beta-distributed cases and controls at a certain significance level
}
\description{mu1_curve_data() calculates power for a range of mu1 values (the mean of the case group), holding mu0 (the mean of the control group) constant, and assuming case and control measurements are beta-distributed. Outputs a dataframe with four columns: mu1, mu0, power_p (the power calculated as a fraction of the p-values below a threshold), power_adj (the power calculated as a fraction of the adjusted p-values below a threshold).
}
\usage{
mu1_curve_data(mu0, sd0, mu1_min, mu1_max, mu1_step, sd1, n0, n1, N, a, MHT, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{a reasonable estimate of the control group's mean
}
\item{sd0}{a reasonable estimate of the control group's standard deviation
}
\item{mu1_min}{a reasonable minimum estimate of the case group's mean
}
\item{mu1_max}{a reasonable maximum estimate of the case group's mean
}
\item{mu1_step}{a reasonable step to increase the mean between the minimum and maximum
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
#Calculate power for a range of effect sizes
data <- mu1_curve_data(mu0=0.2, sd0=0.15, mu1_min=0.25,
mu1_max=0.27, mu1_step=0.01, sd1=0.15,
n0=150, n1=300, N=100, a=0.05,
MHT="bonferroni",q=0.1)
}
