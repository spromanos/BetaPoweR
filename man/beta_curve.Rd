\name{beta_curve}
\alias{beta_curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate power to detect a range of effect sizes between beta-distributed cases and controls at variable matching ratios at a certain significance level
}
\description{beta_curve() calculates the power to detect a range of effect sizes between beta-distributed cases and controls at variable matching ratios at a certain significance level. Outputs a dataframe with six columns: mu1, mu0, power_p (the power calculated as a fraction of the p-values below a threshold), power_adj (the power calculated as a fraction of the adjusted p-values below a threshold), effect_size (the difference of means: mu1-m0), and sample_size (a factor of case:control matching ratios).
}
\usage{
beta_curve(mu0, sd0, mu1_min, mu1_max, mu1_step, sd1, n0s, n1s, N, a, MHT, q)
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
\item{n0s}{a vector of sample sizes of the control group. Each position will be matched to the corresponding position in the n1s vector.
}
\item{n1s}{a vector of sample sizes of the case group. Each position will be matched to the corresponding position in the n0s vector.
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
#Calculate power for a range of effect sizes and variable matching ratios
data <- beta_curve(mu0=0.2, sd0=0.15, mu1_min=0.25, mu1_max=0.27,
mu1_step=0.01, sd1=0.15, n0s=c(150, 150, 150), n1s=c(150, 300, 450),
N=100,a=0.05,MHT="bonferroni",q=0.1)

}
