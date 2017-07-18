### test-ggSurv.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: jul 11 2017 (18:21) 
## Version: 
## last-updated: jul 11 2017 (18:31) 
##           By: Brice
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(lava)

set.seed(10)
n <- 500  

newdata <- data.frame(X1=1)
time <- 0.25

m <- lvm()
regression(m) <- y ~ 1
regression(m) <- s ~ exp(-2*X1)
distribution(m,~X1) <- binomial.lvm()
distribution(m,~cens) <- coxWeibull.lvm(scale=1)
distribution(m,~y) <- a <- coxWeibull.lvm(scale=1,shape=~s)
eventTime(m) <- eventtime ~ min(y=1,cens=0)
d <- as.data.table(sim(m,n))
setkey(d, eventtime)


m.cox <- coxph(Surv(eventtime, status) ~ X1, data = d, y = TRUE, x = TRUE)
ggSurv(m.cox)

mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), data = d, y = TRUE, x = TRUE)
ggSurv(mStrata.cox)

#----------------------------------------------------------------------
### test-ggSurv.R ends here
