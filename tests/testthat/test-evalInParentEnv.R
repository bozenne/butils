### test-evalInParentEnv.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (13:20) 
## Version: 
## last-updated: okt  5 2017 (13:23) 
##           By: Brice Ozenne
##     Update #: 5
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
res1 <- butils:::evalInParentEnv(e.lm$call$data, envir = environment())

df <- data.frame(Y=1:5,X=1:5)
e.lm <- lm(Y~X, data = df)
res2 <- butils:::evalInParentEnv(e.lm$call$data, envir = environment())

test_that("df vs. data.frame(...) ", {
    expect_equal(res1,res2)
})

#----------------------------------------------------------------------
### test-evalInParentEnv.R ends here
