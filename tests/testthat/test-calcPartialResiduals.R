### test-calcPartialResiduals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 25 2017 (20:41) 
## Version: 
## last-updated: maj 26 2017 (15:18) 
##           By: Brice Ozenne
##     Update #: 44
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

context("#### partial residuals ####")

library(lava)
library(nlme)
library(lme4)
library(merTools)
library(AICcmodavg)

# {{{ linear model

set.seed(10)
n <- 100
dd <- data.frame(x0 = rnorm(n),
                 x1 = seq(-3,3, length.out=n),
                 x2 = factor(rep(c(1,2),each=n/2), labels=c("A","B"))
                 )
dd$y <- 5 + 2*dd$x0 + 0.5*dd$x1 + -1*(dd$x2=="B")*dd$x1 + 0.5*(dd$x2=="B") + rnorm(n, sd=0.25)

lm0 <- lm(y ~ x0 + x1*x2, dd)

data(iris) 
l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)

test_that("linear model - 1 variable", {
    res1 <- calcPartialResiduals(lm0, var=c("x1"))
    res2 <- calcPartialResiduals(lm0, var=c("x1"), predictFUN = "predict")
    resGS <- plotConf(lm0, var1 = "x1", plot = FALSE)
    
    expect_equal(res1,res2)

    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    
    expect_equal(as.double(resGS$x),res1$data$x1)
    ## expect_equal(as.double(resGS$y),res1$data$pResiduals)
    ## expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    ## expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    ## expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
    # dev.new()
    # plot(res1)
    # plotConf(lm0,var1="x1")

    res1 <- calcPartialResiduals(l, var=c("Sepal.Width"))
    res2 <- calcPartialResiduals(l, var=c("Sepal.Width"), predictFUN = "predict")
    resGS <- plotConf(l,var1="Sepal.Width", plot = FALSE)
    # dev.new()
    # plot(res1)
    # plotConf(l,var1="Sepal.Width")

    expect_equal(res1,res2)
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = l$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)

    # expect_equal(as.double(resGS$x),res1$data[["Sepal.Width"]])
    ## expect_equal(as.double(resGS$y),res1$data$pResiduals)
    ## expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    ## expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    ## expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
})

test_that("linear model - 2 variables", {
    res1 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    res2 <- calcPartialResiduals(lm0, var=c("x1","x2"), predictFUN = "predict")
    resGS <- plotConf(lm0, var1 = "x1", var2 = "x2", plot = FALSE)

    expect_equal(res1,res2)
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)

    expect_equal(as.double(resGS$x),res1$data$x1)
    ## expect_equal(as.double(resGS$y),res1$data$pResiduals)
    ## expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    ## expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    ## expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)

    res1 <- calcPartialResiduals(l, var=c("Sepal.Width","Species"))
    res2 <- calcPartialResiduals(l, var=c("Sepal.Width","Species"), predictFUN = "predict")
    resGS <- plotConf(l,var1="Sepal.Width",var2="Species", plot = FALSE)
    # dev.new()
    # plot(res1)
    # plotConf(l,var1="Sepal.Width",var2="Species")
    # iris

    expect_equal(res1,res2)
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = l$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    # expect_equal(as.double(resGS$x),res1$data[["Sepal.Width"]])
    ## expect_equal(as.double(resGS$y),res1$data$pResiduals)
    ## expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    ## expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    ## expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
})

# }}}

# {{{ linear model no intercept
set.seed(10)
n <- 100
dd <- data.frame(x0 = rnorm(n),
                 x1 = seq(-3,3, length.out=n),
                 x2 = factor(rep(c(1,2),each=n/2), labels=c("A","B"))
                 )
dd$strata <- as.factor(1:2)
dd$y <- 5 + 2*dd$x0 + 0.5*dd$x1 + -1*(dd$x2=="B")*dd$x1 + 0.5*(dd$x2=="B") + rnorm(n, sd=0.25)
dd.center <- dd
dd.center[,"y"] <- dd[,"y"]-mean(dd[,"y"])

test_that("linear model - 1 variables no intercept", {
    ## only one continuous variable
    lm0 <- lm(y ~ 0 + x0, dd)
    lm0.center <- lm(y ~ 0 +  x0, data = dd.center)
    res1 <- calcPartialResiduals(lm0, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, res1$partialFit))["(Intercept)"]
    expect_equal(as.double(intercept),0)

    res1 <- calcPartialResiduals(lm0.center, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, res1$partialFit))["(Intercept)"]
    expect_equal(as.double(intercept),0)
    # plot(res1)

    ## only one categorical variable
    lm0 <- lm(y ~ 0 + x2, dd)
    res1 <- calcPartialResiduals(lm0, var=c("x2"))
    expect_equal(res1$data$x2, dd$x2)
    expect_equal(res1$data$pResiduals, dd$y)
    # plot(res1)
})

test_that("linear model - 2 variables", {
    ## two categorical
    lm0 <- lm(y ~ 0 + strata + x2, dd)
    res1 <- calcPartialResiduals(lm0, var=c("x2","strata"))
    expect_equal(res1$data$x2, dd$x2)
    expect_equal(res1$data$pResiduals, dd$y)
    # plot(res1)
    
    ## one categorical and one continuous
    res1 <- calcPartialResiduals(lm(y ~ 0 + x0 + x2, dd), var=c("x2"))
    res2 <- calcPartialResiduals(lm(y ~ 0 + x2 + x0, dd), var=c("x2"))
    expect_equal(res1$data$pResiduals,res2$data$pResiduals)
})
# }}}

# {{{ lmer
set.seed(10)
dd$Id <- rbinom(n, size = 3, prob = 0.3)
lm0 <- lm(y ~ x0 + x1*x2, dd)
lmer0 <- lmer(y ~ x0 + x1*x2 + (1|Id), dd)


test_that("lmer", {
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lmer0, var=c("x2"))
    )
    suppressWarnings(
    res2 <- calcPartialResiduals(lmer0, var=c("x2"),
                                 FUN.predict = "predict_merTools")
    )
    suppressWarnings(
    res3 <- calcPartialResiduals(lmer0, var=c("x2"),
                                 FUN.predict = "predict_AICcmodavg")
    )
    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-10)
    #  plot(res0)
    expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-10)
    expect_equal(res1$partialFit$fit, res2$partialFit$fit, tol = 1e-2)
    expect_equal(res1$partialFit$fit, res3$partialFit$fit)
    expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res2$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)

    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lmer0, var=c("x1","x2"))
    )
    suppressWarnings(
    res2 <- calcPartialResiduals(lmer0, var=c("x1","x2"),
                                 FUN.predict = "predict_merTools")
    )
    suppressWarnings(
    res3 <- calcPartialResiduals(lmer0, var=c("x1","x2"),
                                 FUN.predict = "predict_AICcmodavg")
    )
    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-10)
    # plot(res1)
    expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-10)
    expect_equal(res1$partialFit$fit, res2$partialFit$fit, tol = 1e-2)
    expect_equal(res1$partialFit$fit, res3$partialFit$fit)
    expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res2$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)
})

lm0 <- lm(y ~ 0 + x0 + x1*x2, dd)
lmer0 <- lmer(y ~ 0 + x0 + x1*x2 + (1|Id), dd)

test_that("lmer no intercept", {
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lmer0, var=c("x2"))
    )
    #### ERROR
    ## res2 <- calcPartialResiduals(lmer0, var=c("x2"),
    ##                              FUN.predict = "predict_merTools")
    ## res3 <- calcPartialResiduals(lmer0, var=c("x2"),
    ##                              FUN.predict = "predict_AICcmodavg")
    expect_equal(res0$data$pResiduals, res1$data$pResiduals)
    #  plot(res0) 
    
    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    res1 <- calcPartialResiduals(lmer0, var=c("x1","x2"))
    suppressWarnings(
        expect_equal(res0$data$pResiduals, res1$data$pResiduals)
    )
    # plot(res1)
})

# }}}


# {{{ lme
set.seed(10)
dd$Id <- rbinom(n, size = 3, prob = 0.3)
lm0 <- lm(y ~ x0 + x1*x2, dd)
lme0 <- lme(y ~ x0 + x1*x2, random = ~1|Id, dd)

test_that("lme", {
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lme0, var=c("x2"))
    )
    suppressWarnings(
        res3 <- calcPartialResiduals(lme0, var=c("x2"),
                                     FUN.predict = "predict_AICcmodavg")
    )
    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-10)
    #  plot(res0)
    expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-10)
    expect_equal(res1$partialFit$fit, res3$partialFit$fit)
    expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)

    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lme0, var=c("x1","x2"))
    )
    suppressWarnings(
        res3 <- calcPartialResiduals(lme0, var=c("x1","x2"),
                                     FUN.predict = "predict_AICcmodavg")
    )
    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-10)
    # plot(res1)
    expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-10)
    expect_equal(res1$partialFit$fit, res3$partialFit$fit)
    expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
    expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)
})

lm0 <- lm(y ~ 0 + x0 + x1*x2, dd)
lme0 <- lme(y ~ 0 + x0 + x1*x2, random = ~ 1|Id, dd)

test_that("lme no intercept", {
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    suppressWarnings(
        res1 <- calcPartialResiduals(lme0, var=c("x2"))
    )
    expect_equal(res0$data$pResiduals, res1$data$pResiduals)
    #  plot(res0) 
    
    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    res1 <- calcPartialResiduals(lme0, var=c("x1","x2"))
    suppressWarnings(
        expect_equal(res0$data$pResiduals, res1$data$pResiduals)
    )
    # plot(res1)
})

# }}}


## # {{{ gls
## set.seed(10)
## n <- 100
## dd <- data.frame(x0 = rnorm(n),
##                  x1 = seq(-3,3, length.out=n),
##                  x2 = factor(rep(c(1,2),each=n/2), labels=c("A","B"))
##                  )
## dd$y <- 5 + 2*dd$x0 + 0.5*dd$x1 + -1*(dd$x2=="B")*dd$x1 + 0.5*(dd$x2=="B") + rnorm(n, sd=0.25)
## dd$Id <- rbinom(n, size = 3, prob = 0.3)
## lm0 <- lm(y ~ x0 + x1*x2, dd)
## gls0 <- gls(y ~ x0 + x1*x2, correlation = corCompSymm(form=~ 1|Id), dd)

## test_that("gls", {
##     res0 <- calcPartialResiduals(lm0, var=c("x2"))
##     suppressWarnings(
##         res1 <- calcPartialResiduals(gls0, var=c("x2"))
##     )
##     suppressWarnings(
##         res3 <- calcPartialResiduals(gls0, var=c("x2"),
##                                      FUN.predict = "predict_AICcmodavg")
##     )
##     expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-2)
##     #  plot(res0)
##     expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-2)
##     expect_equal(res1$partialFit$fit, res3$partialFit$fit)
##     expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
##     expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)

##     res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
##     suppressWarnings(
##         res1 <- calcPartialResiduals(gls0, var=c("x1","x2"))
##     )
##     suppressWarnings(
##         res3 <- calcPartialResiduals(gls0, var=c("x1","x2"),
##                                      FUN.predict = "predict_AICcmodavg")
##     )
##     expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-3)
##     # plot(res1)
##     expect_equal(res1$partialFit$fit, res0$partialFit$fit, tol = 1e-2)
##     expect_equal(res1$partialFit$fit, res3$partialFit$fit)
##     expect_equal(res1$partialFit$fit.lower, res0$partialFit$fit.lower, tol = 1e-1)
##     expect_equal(res1$partialFit$fit.lower, res3$partialFit$fit.lower)
## })

## lm0 <- lm(y ~ 0 + x0 + x1*x2, dd)
## gls0 <- gls(y ~ 0 + x0 + x1*x2, correlation = corCompSymm(form=~ 1|Id), dd)

## test_that("gls no intercept", {
##     res0 <- calcPartialResiduals(lm0, var=c("x2"))
##     suppressWarnings(
##         res1 <- calcPartialResiduals(gls0, var=c("x2"))
##     )
##     expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-3)
##     #  plot(res0) 
    
##     res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
##     res1 <- calcPartialResiduals(gls0, var=c("x1","x2"))
##     suppressWarnings(
##         expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-3)
##     )
##     # plot(res1)
## })

## # }}}

#----------------------------------------------------------------------
### test-calcPartialResiduals.R ends here
