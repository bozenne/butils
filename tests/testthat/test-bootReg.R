# library(butils.base)
# package.source("butils")
library(testthat)

context("#### bootReg ####")

## * linear model
library(lava)

set.seed(10)
n <- 1e2
mSim <- lvm()
regression(mSim, Y~X1+X2+X3+X4+X5) <- c(2,1,0,0,-1)
df.data <- lava::sim(mSim, n)

e.lm <- lm(Y~X1+X2+X3+X4+X5, data = df.data)
boot.lm <- bootReg(e.lm, n.boot = 5e2)

## compare
sboot.lm <- summary(boot.lm)
sboot.lm
summary(e.lm)$coef


## * under H1

test.sim <- FALSE

if(test.sim){

    library(lmeresampler)
  
    vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class,
                  random = ~1 | school, data = jsp728)
  
    ## you can write your own function to return stats, or use something like 'fixef'
    mySumm <- function(.) { 
        return(fixef(., "beta")) 
    }

    
    ## running a parametric bootstrap 
    n.boot <- 1e3
    system.time(
        resPackage <- bootstrap(model = vcmodA, fn = mySumm, type = "parametric", B = n.boot)
    )
    system.time(
        resButils <- bootReg(vcmodA, n.boot = n.boot)
    )
    s.resButils <- summary(resButils)
    s.resButils
    summary(vcmodA)$tTable
    
    fm1 <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
               random = ~ 1 | Mare)
  
}

