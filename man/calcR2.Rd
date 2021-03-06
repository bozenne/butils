% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R2.R
\name{calcR2}
\alias{calcR2}
\title{Compute the R square using the predicted values}
\usage{
calcR2(model, data = NULL, trace = FALSE)
}
\arguments{
\item{model}{the model from which the R squared should be computed.}

\item{data}{the data that have been used to fit the model.}

\item{trace}{should the execution of the function be traced.}
}
\description{
Compute the R square using the predicted values
}
\details{
Compute the Generalized R2 defined by Buse (1973) using the formula:
\deqn{R^2 = 1- \frac{e S^{-1} e}{e0 S^{-1} e0}}
where iS is the inverse of the variance-covariance matrix of the observations.
      e is the vector of residuals of the full model.
      e0 is the vector of residuals of the reduced model.
      
Compute the Generalized R2 based on the log-likelihood:
\deqn{R^2 = 1 - exp(- LR / n)}
where n is the number of observations.
      LR is twice the log of the ratio of the likelihood between the full and reduced model.

Denoting X the covariate of interest, Y the outcome, and Z the other covariates, 
the partial R^2 is the same as:
\itemize{
\item the square of the partial correlation coefficient between X and Y, controlling for Z.
\item the partial \eqn{\eta^2} which is the ratio between the sum of squares from X (SSX) and the sum of SSX and the residual sum of squares (SSE)
}
}
\examples{
library(lava)
library(nlme)

m <- lvm(Y~X1+X2+X3, G~1)
categorical(m, K=5, label = c("A","B","C","D","E")) <- ~G
d <- lava::sim(m, 1e2)
d <- d[order(d$G),]
d$Y <- d$Y + 0.5*as.numeric(as.factor(d$G))

## linear model
m.lm <- lm(Y~X1+X2+X3, data = d)
calcR2(m.lm)
summary(m.lm)

summary(m.lm)$r.squared

ff <- Y~X1+X2+X3
m.lm <- lm(ff, data = d)
calcR2(m.lm)
summary(m.lm)

dt <- as.data.table(d)
dt[,interaction := X1*X2]
m.lm2 <- lm(Y~X1+X2+X3+interaction, data = dt[G \%in\% c("A","B")])
calcR2(m.lm2)

## gls model
m.gls1 <- gls(Y~X1+X2+X3, data = d)
calcR2(m.gls1)

## heteroschedasticity
m.gls2 <- gls(Y~X1+X2+X3,
              weights = varIdent(form = ~ 1 |G), data = d)
calcR2(m.gls2)

## correlation
m.gls2 <- gls(Y~X1+X2+X3,
              correlation = corCompSymm(form = ~ 1 |G), data = d)
calcR2(m.gls2)

}
\references{
A. Buse (1973). Goodness of Fit in Generalized Least Squares Estimation. The American Statistician,, Vol. 27, No. 3.
Lonnie Magee (1990). R2 Measures Based on Wald and Likelihood Ratio Joint Significance Tests, The American Statistician, 44:3, 250-253.
https://stats.stackexchange.com/questions/64010/importance-of-predictors-in-multiple-regression-partial-r2-vs-standardized
}
