% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditionalCI.R
\name{conditionalCI}
\alias{conditionalCI}
\title{Conditional confidence interval}
\usage{
conditionalCI(
  theta,
  threshold,
  sigma = 1,
  conf.level = 0.95,
  df = NULL,
  method = "shortest",
  distribution = "gaussian",
  trace = length(theta) > 1,
  ...
)
}
\arguments{
\item{theta}{[numeric vector] the observed value(s) of the statistic.}

\item{threshold}{[numeric, >0] the threshold.}

\item{sigma}{[numeric vector] the standard error of the statistic(s).}

\item{conf.level}{[numeric, 0-1] Conditional confidence level of the interval.}

\item{df}{[interger, >0] Degree of freedom of the distribution of theta.
Only relevant when distribution equals \code{"student"}.}

\item{method}{[character] the method used to compute the conditional confidence interval.
So far only \code{"shortest"} is supported.}

\item{distribution}{[character] Distribution of theta.
Either \code{"gaussian"} or \code{"student"}.}

\item{trace}{[interger, 0-2] should a progress bar be displayed?}

\item{...}{additional arguments to be passed to \code{lavaSearch2:::.calcShortestCI}
to specify the optimization method.}
}
\description{
Conditional confidence intervals for normally distributed random variables.
}
\details{
Compute the confidence interval of theta
conditional on theta being greater (in absolute value) than the threshold.
theta is assumed to be normally distributed with known variance sigma.
The thresholding is applied after normalization (i.e. dividing theta by sigma).
This corresponding to conditioning on a minimal significance level.
}
\examples{
ci <- conditionalCI(theta = seq(3,7, by = 0.1),
                    threshold = 2.9999, distribution = "gaussian")

print(ci)
confint(ci)
autoplot(ci)

ci1 <- conditionalCI(theta = 4, sigma = 2, threshold = 2.9999, method = "shortest")
ci2 <- conditionalCI(theta = 4/2, sigma = 1, threshold = 2.9999/2, method = "shortest")
confint(ci1) - 2*confint(ci2)
ci3 <- conditionalCI(theta = 4, sigma = 2, threshold = 2.9999, method = "shortest2")
confint(ci1) - confint(ci3)
}
