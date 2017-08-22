### test-ggSurv.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: jul 11 2017 (18:21) 
## Version: 
## last-updated: jul 18 2017 (11:55) 
##           By: Brice Ozenne
##     Update #: 11
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


context("#### Display survival curves ####")

# {{{ random example
library(survival)
dt <- as.data.table(aml)
dt[,x:=as.factor(x)]
 
## survfit
KM <- survfit(Surv(time, status) ~ x, data = dt)
ggSurv(KM)
ggSurv(KM, censoring = TRUE, event = TRUE)
ggSurv(KM, confint = TRUE)
  
## survfit
KM <- survfit(Surv(time, status) ~ x, data = dt)
ggSurv(KM)
ggSurv(KM, censoring = TRUE, event = TRUE)
ggSurv(KM, confint = TRUE)
 
## data.table
dt2 <- data.table(time = 1:10, 
                  survival = seq(1, by = -0.01, length.out = 10), 
                  n.censor = 0, n.event = 1)
ggSurv(dt2)

dt3 <- data.table(time = 1:10, 
                  survival = seq(1, by = -0.01, length.out = 10)
                  )
ggSurv(dt3)
# }}} 

# {{{ display of events and censoring
d <- data.table(time = 1:4,
                status = rep(0:1,2)
                )
KM <- survfit(Surv(time, status) ~ 1, data = d)
ggSurv(KM, censoring = TRUE, event = TRUE)

dS <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2),
            cbind(time = 0.5, status = 1, strata = 2))
KM <- survfit(Surv(time, status) ~ strata, data = dS)
ggSurv(KM, censoring = TRUE, event = TRUE)

# }}}
 
#----------------------------------------------------------------------
### test-ggSurv.R ends here
