### test-calcPartialResiduals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 25 2017 (20:41) 
## Version: 
## last-updated: jan 24 2018 (17:31) 
##           By: Brice Ozenne
##     Update #: 70
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

## * dataset
set.seed(10)
n <- 100
dd <- data.frame(x0 = rnorm(n),
                 x1 = seq(-3,3, length.out=n),
                 x2 = factor(rep(c(1,2),each=n/2), labels=c("A","B"))
                 )
dd$strata <- as.factor(1:2)
dd$Id <- rbinom(n, size = 3, prob = 0.3)
dd$y <- 5 + 2*dd$x0 + 0.5*dd$x1 + -1*(dd$x2=="B")*dd$x1 + 0.5*(dd$x2=="B") + rnorm(n, sd=0.25)
dd$y.center <- dd$y-mean(dd$y)

## * linear model
test_that("linear model (1) - 1 variable", {
    lm0 <- lm(y ~ x0 + x1*x2, data = dd)
  
    res1 <- calcPartialResiduals(lm0, var=c("x1"), keep.intercept = TRUE)
    res2 <- calcPartialResiduals(lm0, var=c("x1"), keep.intercept = FALSE)
    resGS <- plotConf(lm0, var1 = "x1", plot = FALSE)
    
    expect_equal(res1$data$pResiduals, as.double(resGS$y))
    expect_equal(res1$data$pResiduals,  res2$data$pResiduals + coef(lm0)[1])
    
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    
    expect_equal(as.double(resGS$x),res1$partialFit$x1)
    expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
    expect_equal(res1$partialFit$fit,res2$partialFit$fit + coef(lm0)[1])
})

test_that("linear model (2) - 1 variable", {
    data(iris) 
    lm0 <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
    
    res1 <- calcPartialResiduals(lm0, var=c("Sepal.Width"), keep.intercept = TRUE)
    res2 <- calcPartialResiduals(lm0, var=c("Sepal.Width"), keep.intercept = FALSE)
    resGS <- plotConf(lm0,var1="Sepal.Width", plot = FALSE)
    
    expect_equal(res1$data$pResiduals, as.double(resGS$y))
    expect_equal(res1$data$pResiduals,  res2$data$pResiduals + coef(lm0)[1])
    
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)

    expect_equal(as.double(resGS$x),res1$partialFit$Sepal.Width)
    expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
    expect_equal(res1$partialFit$fit,res2$partialFit$fit + coef(lm0)[1])
})

test_that("linear model (1) - 2 variables", {
    lm0 <- lm(y ~ x0 + x1*x2, data = dd)

    res1 <- calcPartialResiduals(lm0, var=c("x1","x2"), keep.intercept = TRUE)
    res2 <- calcPartialResiduals(lm0, var=c("x1","x2"), keep.intercept = FALSE)
    resGS <- plotConf(lm0, var1="x1", var2 ="x2", plot = FALSE)
    
    expect_equal(res1$data$pResiduals, as.double(resGS$y))
    expect_equal(res1$data$pResiduals,  res2$data$pResiduals + coef(lm0)[1])
    
    
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)

    expect_equal(as.double(rep(resGS$x,2)),res1$partialFit$x1)
    expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
    expect_equal(res1$partialFit$fit,res2$partialFit$fit + coef(lm0)[1])
})

test_that("linear model (2) - 2 variables", {
    data(iris) 
    lm0 <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
    
    res1 <- calcPartialResiduals(lm0, var=c("Sepal.Width","Species"), keep.intercept = TRUE)
    res2 <- calcPartialResiduals(lm0, var=c("Sepal.Width","Species"), keep.intercept = FALSE)
    resGS <- plotConf(lm0, var1="Sepal.Width", var2 ="Species", plot = FALSE)
    
    expect_equal(res1$data$pResiduals, as.double(resGS$y))
    expect_equal(res1$data$pResiduals,  res2$data$pResiduals + coef(lm0)[1])
    
    
    # to get normal quantile instead of student to match lava
    cst <- qnorm(1 - (1 - 0.95)/2)/qt(1 - (1 - 0.95)/2, df = lm0$df.residual)
    res1$partialFit$fit.lower <- res1$partialFit$fit - cst*(res1$partialFit$fit.upper-res1$partialFit$fit)
    res1$partialFit$fit.upper <- res1$partialFit$fit + cst*(res1$partialFit$fit.upper-res1$partialFit$fit)

    expect_equal(as.double(rep(resGS$x,3)),res1$partialFit$Sepal.Width)
    expect_equal(as.double(resGS$predict$fit[,1]),res1$partialFit$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),res1$partialFit$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),res1$partialFit$fit.upper)
    expect_equal(res1$partialFit$fit,res2$partialFit$fit + coef(lm0)[1])

})

## * linear model no intercept
test_that("linear model - 1 variables no intercept", {
    ## only one continuous variable
    lm0 <- lm(y ~ 0 + x0, data = dd)
    lm0.center <- lm(y.center ~ 0 +  x0, data = dd)
    
    res1 <- calcPartialResiduals(lm0, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, res1$partialFit))["(Intercept)"]
    expect_equal(as.double(intercept),0)

    res1 <- calcPartialResiduals(lm0.center, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, res1$partialFit))["(Intercept)"]
    expect_equal(as.double(intercept),0)
    # plot(res1)

    ## only one categorical variable
    lm0 <- lm(y ~ 0 + x2, data = dd)
    res1 <- calcPartialResiduals(lm0, var=c("x2"), keep.intercept = TRUE)
    expect_equal(res1$data$x2, dd$x2)
    expect_equal(res1$data$pResiduals, dd$y)
    # plot(res1)
})

test_that("linear model - 2 variables", {
    ## two categorical
    lm0 <- lm(y ~ 0 + strata + x2, dd)
    res1 <- calcPartialResiduals(lm0, var=c("x2","strata"), keep.intercept = TRUE)
    expect_equal(res1$data$x2, dd$x2)
    expect_equal(res1$data$pResiduals, dd$y)
    # plot(res1)
    
    ## one categorical and one continuous
    res1 <- calcPartialResiduals(lm(y ~ 0 + x0 + x2, dd), var=c("x2"))
    res2 <- calcPartialResiduals(lm(y ~ 0 + x2 + x0, dd), var=c("x2"))
    expect_equal(res1$data$pResiduals,res2$data$pResiduals)
})


## * mixed models

## ** no random effects
test_that("no random effect", {
    lm0 <- lm(y ~ x0 + x1*x2, dd)
    lmer0 <- lmer(y ~ x0 + x1*x2 + (1|Id), dd)
    lme0 <- lme(y ~ x0 + x1*x2, random = ~1|Id, dd)
    
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    res1 <- calcPartialResiduals(lmer0, var=c("x2"))
    res2 <- calcPartialResiduals(lme0, var=c("x2"))

    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res1$partialFit$fit, tol = 1e-6)
    expect_equal(res0$data$pResiduals, res2$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res2$partialFit$fit, tol = 1e-6)

    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    res1 <- calcPartialResiduals(lmer0, var=c("x1","x2"))
    res2 <- calcPartialResiduals(lme0, var=c("x1","x2"))

    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res1$partialFit$fit, tol = 1e-6)
    expect_equal(res0$data$pResiduals, res2$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res2$partialFit$fit, tol = 1e-6)
})


test_that("no random effect no intercept", {
    lm0 <- lm(y ~ 0 + x0 + x1*x2, dd)
    lmer0 <- lmer(y ~ 0 + x0 + x1*x2 + (1|Id), dd)
    lme0 <- lme(y ~ 0 + x0 + x1*x2, random = ~1|Id, dd)
    
    res0 <- calcPartialResiduals(lm0, var=c("x2"))
    res1 <- calcPartialResiduals(lmer0, var=c("x2"))
    res2 <- calcPartialResiduals(lme0, var=c("x2"))

    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res1$partialFit$fit, tol = 1e-6)
    expect_equal(res0$data$pResiduals, res2$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res2$partialFit$fit, tol = 1e-6)

    res0 <- calcPartialResiduals(lm0, var=c("x1","x2"))
    res1 <- calcPartialResiduals(lmer0, var=c("x1","x2"))
    res2 <- calcPartialResiduals(lme0, var=c("x1","x2"))

    expect_equal(res0$data$pResiduals, res1$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res1$partialFit$fit, tol = 1e-6)
    expect_equal(res0$data$pResiduals, res2$data$pResiduals, tol = 1e-6)
    expect_equal(res0$partialFit$fit, res2$partialFit$fit, tol = 1e-6)
})


## *  gls
test_that("gls", {
    gls0 <- gls(y ~ x0 + x1*x2, correlation = corCompSymm(form=~ 1|Id), dd)

    gls0.X <- model.matrix(y ~ x0 + x1*x2,dd)
    
    res1 <- calcPartialResiduals(gls0,  var=c("x0"), keep.intercept = FALSE)
    res2 <- calcPartialResiduals(gls0, var=c("x0"), keep.intercept = TRUE)
    GS.F <- as.double(gls0.X[,"x0"] * coef(gls0)["x0"] + residuals(gls0))
    GS.T <- as.double(coef(gls0)["(Intercept)"] +  gls0.X[,"x0"] * coef(gls0)["x0"] + residuals(gls0))

    expect_equal(GS.F, res1$data$pResiduals)
    expect_equal(GS.T, res2$data$pResiduals)    
})


#----------------------------------------------------------------------
### test-calcPartialResiduals.R ends here
