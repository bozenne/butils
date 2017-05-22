# library(butils.base)
# package.source("butils")
library(testthat)

context("#### bootGLS ####")

# {{{ linear model

set.seed(10)
n <- 1e2
n.boot <- 10
n.sim <- 2
diff <- 0.3


ls.pvalue <- list(H0=NULL,H1=NULL)
ls.CI <- list(H0=NULL,H1=NULL)
pb <- utils:::txtProgressBar(max =  n.sim)
for(iSim in 1:n.sim){ # iSim <- 1

    for(iDiff in 1:2){ # iDiff <- 1
        Y1 <- rnorm(n, mean = 0)
        Y2 <- rnorm(n, mean = c(0,diff)[iDiff])
        df <- rbind(data.frame(Y=Y1,G=1),
                    data.frame(Y=Y2,G=2)            
                    )
        m.lm <- lm(Y ~ G, data = df)

        resH1 <- bootGLS(m.lm, n.boot = n.boot)
        resH0 <- bootGLS(m.lm, n.boot = n.boot, GROUPvar = "G")

        ls.CI[[iDiff]] <- rbind(ls.CI[[iDiff]],
                                c(resH0$p.value["G"]>0.05,resH1$p.value["G"]>0.05)
                                )
        ls.pvalue[[iDiff]] <- rbind(ls.pvalue[[iDiff]],
                                    c( (resH0$CI["2.5%","G"]<=diff)*(resH0$CI["97.5%","G"]>=diff),
                                    (resH1$CI["2.5%","G"]<=diff)*(resH1$CI["97.5%","G"]>=diff)
                                    )
                                    )
    }
    utils:::setTxtProgressBar(pb, iSim)
}
close(pb)

# }}}

# {{{ under H1

# }}}

# {{{
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
    resButils <- bootGLS(vcmodA, n.boot = n.boot)
  )
  
  
  (apply(resButils$all.boot,2,quantile)-apply(resPackage$t,2,quantile))/apply(resPackage$t,2,quantile)
  
  
  fm1 <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
              random = ~ 1 | Mare)
  
  ## running a parametric bootstrap 
  n.boot <- 1e3
  system.time(
    resPackage <- bootstrap(model = fm1, fn = mySumm, type = "parametric", B = n.boot)
  )
  system.time(
    resButils <- bootGLS(fm1, n.boot = n.boot)
  )
  
}
# }}}