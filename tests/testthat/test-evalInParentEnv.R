### test-evalInParentEnv.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (13:20) 
## Version: 
## last-updated: sep 28 2018 (13:45) 
##           By: Brice Ozenne
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
context("evalInParentEnv")


e.lm <- lm(Y~X, data = data.frame(Y=1:5,X=1:5))
res1 <- butils:::evalInParentEnv(e.lm$call$data)

df <- data.frame(Y=1:5,X=1:5)
e.lm <- lm(Y~X, data = df)
res2 <- butils:::evalInParentEnv(e.lm$call$data)

test_that("df vs. data.frame(...) ", {
    expect_equal(res1,res2)
})

#----------------------------------------------------------------------
### test-evalInParentEnv.R ends here
