% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoplot.partialResiduals.R
\name{autoplot.partialResiduals}
\alias{autoplot.partialResiduals}
\title{Display an first order interaction for categorical variables}
\usage{
\method{autoplot}{partialResiduals}(object, alpha.ggplot = 0.25, plot = TRUE, ...)
}
\arguments{
\item{alpha.ggplot}{transparency parameter for the points or the bands}

\item{plot}{should the plot be displayed?}

\item{...}{additional argument to be passed to \code{partialModel}}

\item{x}{a linear model}
}
\description{
Display an first order interaction for categorical variables
}
\examples{
library(lava)
set.seed(10)
m.lvm <- lvm(Y~2*Age+4*gender+gene+time)
categorical(m.lvm, labels = c("M","F")) <- ~gender
categorical(m.lvm, K = 10) <- ~gene
d <- lava::sim(n = 1e2, m.lvm)

## linear model
m <- lm(Y~Age+gender, data = d)
pres1 <- partialResiduals(m, var = "Age")
autoplot(pres1)
pres2 <- partialResiduals(m, var = c("Age","gender"))
autoplot(pres2)
pres3 <- partialResiduals(m, var = c("Age","gender"), interval = "prediction")
autoplot(pres3)
pres4 <- partialResiduals(m, var = "gender")
autoplot(pres4)

## linear mixed model
if(require(nlme)){
mm <- lme(Y~Age+gender, random= ~ 1|gene, data = d)
pres1 <- partialResiduals(mm, var = "Age")
autoplot(pres1)
}
if(require(lme4) && require(merTools) && require(AICcmodavg)){

mm <- lmer(Y~Age+gender+(1|gene), data = d)
pres1 <- partialResiduals(mm, var = "Age")
autoplot(pres1)
pres2 <- partialResiduals(mm, var = c("Age","gender"))
autoplot(pres2)

# using external function
pres3 <- partialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_merTools)
autoplot(pres3) 
pres4 <- partialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_AICcmodavg)
autoplot(pres4)
}

## gam
if(require(mgcv)){
set.seed(2) ## simulate some data
dat <- gamSim(1,n=400,dist="normal",scale=2)
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
b <- gam(y~s(x0)+x1+x2+x3,data=dat)
pres5 <- partialResiduals(b, var = "x0")
autoplot(pres5)
}
}
