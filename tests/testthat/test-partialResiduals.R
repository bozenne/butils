### test-partialResiduals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 25 2017 (20:41) 
## Version: 
## last-updated: jun 11 2020 (11:21) 
##           By: Brice Ozenne
##     Update #: 78
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

    ## no df to match lava results
    res1 <- partialResiduals(lm0, var=c("x1"), keep.intercept = TRUE, FUN.df = function(...){NULL})
    res2 <- partialResiduals(lm0, var=c("x1"), keep.intercept = FALSE, FUN.df = function(...){NULL})
    resGS <- plotConf(lm0, var1 = "x1", plot = FALSE)
    
    expect_equal(res1$pResiduals, as.double(resGS$y))
    expect_equal(res1$pResiduals,  res2$pResiduals + coef(lm0)[1])
    
    expect_equal(as.double(resGS$x),attr(res1,"partialFit")$x1)
    expect_equal(as.double(resGS$predict$fit[,1]),attr(res1,"partialFit")$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),attr(res1,"partialFit")$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),attr(res1,"partialFit")$fit.upper)
    expect_equal(attr(res1,"partialFit")$fit,attr(res2,"partialFit")$fit + coef(lm0)[1])
})

test_that("linear model (2) - 1 variable", {
    data(iris) 
    lm0 <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
    
    ## no df to match lava results
    res1 <- partialResiduals(lm0, var=c("Sepal.Width"), keep.intercept = TRUE, FUN.df = function(...){NULL})
    res2 <- partialResiduals(lm0, var=c("Sepal.Width"), keep.intercept = FALSE, FUN.df = function(...){NULL})
    resGS <- plotConf(lm0,var1="Sepal.Width", plot = FALSE)
    
    expect_equal(res1$pResiduals, as.double(resGS$y))
    expect_equal(res1$pResiduals,  res2$pResiduals + coef(lm0)[1])
    
    expect_equal(as.double(resGS$x),attr(res1,"partialFit")$Sepal.Width)
    expect_equal(as.double(resGS$predict$fit[,1]),attr(res1,"partialFit")$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),attr(res1,"partialFit")$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),attr(res1,"partialFit")$fit.upper)
    expect_equal(attr(res1,"partialFit")$fit,attr(res2,"partialFit")$fit + coef(lm0)[1])
})

test_that("linear model (1) - 2 variables", {
    lm0 <- lm(y ~ x0 + x1*x2, data = dd)

    ## no df to match lava results
    res1 <- partialResiduals(lm0, var=c("x1","x2"), keep.intercept = TRUE, FUN.df = function(...){NULL})
    res2 <- partialResiduals(lm0, var=c("x1","x2"), keep.intercept = FALSE, FUN.df = function(...){NULL})
    resGS <- plotConf(lm0, var1="x1", var2 ="x2", plot = FALSE)
    
    expect_equal(res1$pResiduals, as.double(resGS$y))
    expect_equal(res1$pResiduals,  res2$pResiduals + coef(lm0)[1])
    
    expect_equal(as.double(rep(resGS$x,2)),attr(res1,"partialFit")$x1)
    expect_equal(as.double(resGS$predict$fit[,1]),attr(res1,"partialFit")$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),attr(res1,"partialFit")$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),attr(res1,"partialFit")$fit.upper)
    expect_equal(attr(res1,"partialFit")$fit,attr(res2,"partialFit")$fit + coef(lm0)[1])
})

test_that("linear model (2) - 2 variables", {
    data(iris) 
    lm0 <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
    
    ## no df to match lava results
    res1 <- partialResiduals(lm0, var=c("Sepal.Width","Species"), keep.intercept = TRUE, FUN.df = function(...){NULL})
    res2 <- partialResiduals(lm0, var=c("Sepal.Width","Species"), keep.intercept = FALSE, FUN.df = function(...){NULL})
    resGS <- plotConf(lm0, var1="Sepal.Width", var2 ="Species", plot = FALSE)
    
    expect_equal(res1$pResiduals, as.double(resGS$y))
    expect_equal(res1$pResiduals,  res2$pResiduals + coef(lm0)[1])
    
    
    expect_equal(as.double(rep(resGS$x,3)),attr(res1,"partialFit")$Sepal.Width)
    expect_equal(as.double(resGS$predict$fit[,1]),attr(res1,"partialFit")$fit)
    expect_equal(as.double(resGS$predict$fit[,2]),attr(res1,"partialFit")$fit.lower)
    expect_equal(as.double(resGS$predict$fit[,3]),attr(res1,"partialFit")$fit.upper)
    expect_equal(attr(res1,"partialFit")$fit,attr(res2,"partialFit")$fit + coef(lm0)[1])

})

## * linear model no intercept
test_that("linear model - 1 variables no intercept", {
    ## only one continuous variable
    lm0 <- lm(y ~ 0 + x0, data = dd)
    lm0.center <- lm(y.center ~ 0 +  x0, data = dd)
    
    res1 <- partialResiduals(lm0, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, attr(res1,"partialFit")))["(Intercept)"]
    expect_equal(as.double(intercept),0)

    res1 <- partialResiduals(lm0.center, var=c("x0"))
    intercept <- coef(lm(fit ~ x0, attr(res1,"partialFit")))["(Intercept)"]
    expect_equal(as.double(intercept),0)
    # plot(res1)

    ## only one categorical variable
    lm0 <- lm(y ~ 0 + x2, data = dd)
    res1 <- partialResiduals(lm0, var=c("x2"), keep.intercept = TRUE)
    expect_equal(res1$x2, dd$x2)
    expect_equal(res1$pResiduals, dd$y)
    # plot(res1)
})

test_that("linear model - 2 variables", {
    ## two categorical
    lm0 <- lm(y ~ 0 + strata + x2, dd)
    res1 <- partialResiduals(lm0, var=c("x2","strata"), keep.intercept = TRUE)
    expect_equal(res1$x2, dd$x2)
    expect_equal(res1$pResiduals, dd$y)
    # plot(res1)
    
    ## one categorical and one continuous
    res1 <- partialResiduals(lm(y ~ 0 + x0 + x2, dd), var=c("x2"))
    res2 <- partialResiduals(lm(y ~ 0 + x2 + x0, dd), var=c("x2"))
    expect_equal(res1$pResiduals,res2$pResiduals)
})


## * mixed models

## ** no random effects
test_that("no random effect", {
    lm0 <- lm(y ~ x0 + x1*x2, dd)
    lmer0 <- lmer(y ~ x0 + x1*x2 + (1|Id), dd)
    lme0 <- lme(y ~ x0 + x1*x2, random = ~1|Id, dd)
    
    res0 <- partialResiduals(lm0, var=c("x2"))
    res1 <- partialResiduals(lmer0, var=c("x2"))
    res2 <- partialResiduals(lme0, var=c("x2"))

    expect_equal(res0$pResiduals, res1$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res1,"partialFit")$fit, tol = 1e-6)
    expect_equal(res0$pResiduals, res2$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res2,"partialFit")$fit, tol = 1e-6)

    res0 <- partialResiduals(lm0, var=c("x1","x2"))
    res1 <- partialResiduals(lmer0, var=c("x1","x2"))
    res2 <- partialResiduals(lme0, var=c("x1","x2"))

    expect_equal(res0$pResiduals, res1$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res1,"partialFit")$fit, tol = 1e-6)
    expect_equal(res0$pResiduals, res2$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res2,"partialFit")$fit, tol = 1e-6)
})


test_that("no random effect no intercept", {
    lm0 <- lm(y ~ 0 + x0 + x1*x2, dd)
    lmer0 <- lmer(y ~ 0 + x0 + x1*x2 + (1|Id), dd)
    lme0 <- lme(y ~ 0 + x0 + x1*x2, random = ~1|Id, dd)
    
    res0 <- partialResiduals(lm0, var=c("x2"))
    res1 <- partialResiduals(lmer0, var=c("x2"))
    res2 <- partialResiduals(lme0, var=c("x2"))

    expect_equal(res0$pResiduals, res1$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res1,"partialFit")$fit, tol = 1e-6)
    expect_equal(res0$pResiduals, res2$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res2,"partialFit")$fit, tol = 1e-6)

    res0 <- partialResiduals(lm0, var=c("x1","x2"))
    res1 <- partialResiduals(lmer0, var=c("x1","x2"))
    res2 <- partialResiduals(lme0, var=c("x1","x2"))

    expect_equal(res0$pResiduals, res1$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res1,"partialFit")$fit, tol = 1e-6)
    expect_equal(res0$pResiduals, res2$pResiduals, tol = 1e-6)
    expect_equal(attr(res0,"partialFit")$fit, attr(res2,"partialFit")$fit, tol = 1e-6)
})


## * gls
test_that("gls", {
    gls0 <- gls(y ~ x0 + x1*x2, correlation = corCompSymm(form=~ 1|Id), dd)

    gls0.X <- model.matrix(y ~ x0 + x1*x2,dd)
    
    res1 <- partialResiduals(gls0,  var=c("x0"), keep.intercept = FALSE)
    res2 <- partialResiduals(gls0, var=c("x0"), keep.intercept = TRUE)
    GS.F <- as.double(gls0.X[,"x0"] * coef(gls0)["x0"] + residuals(gls0))
    GS.T <- as.double(coef(gls0)["(Intercept)"] +  gls0.X[,"x0"] * coef(gls0)["x0"] + residuals(gls0))

    expect_equal(GS.F, res1$pResiduals)
    expect_equal(GS.T, res2$pResiduals)    
})

## * categorical variables and missing values
mSim <- lvm(AUC ~ OC + Age + Work)
categorical(mSim, labels = c("Yes","No")) <- ~OC
categorical(mSim, labels = c("Study","Rest","Work")) <- ~Work

set.seed(10)
d <- lava::sim(mSim, n=1e2)
d[1,] <- NA

test_that("partial residuals (categorical variables and missing values)", {
    ## no interaction
    e <- lm(AUC ~ OC + Age + Work, data = d)
    e.PR <- partialResiduals(e, var = "OC")

    X <- model.matrix(e)
    keep.coef <- c("(Intercept)","Age","WorkRest","WorkWork")
    GS <- na.omit(d)$AUC - X[,keep.coef] %*% coef(e)[keep.coef]
    expect_equal(as.double(GS[,1]), as.double(e.PR$pResiduals))

    ## interaction
    e.I <- lm(AUC ~ OC * Age + Work, data = d)

    e.I.PR <- partialResiduals(e.I, var = "OC")
    X.I <- model.matrix(e.I)
    keep.coef.I <- c("(Intercept)","Age","WorkRest","WorkWork")
    GS.I <- na.omit(d)$AUC - X.I[,keep.coef.I] %*% coef(e.I)[keep.coef.I]
    expect_equal(as.double(GS.I[,1]), as.double(e.I.PR$pResiduals))

    e.I.PR <- partialResiduals(e.I, var = c("OC","Work"))
    keep.coef.I <- c("(Intercept)","Age")
    GS.I <- na.omit(d)$AUC - X.I[,keep.coef.I] %*% coef(e.I)[keep.coef.I]
    expect_equal(as.double(GS.I[,1]), as.double(e.I.PR$pResiduals))
})
#----------------------------------------------------------------------
### test-partialResiduals.R ends here
