# library(butils.base)
# package.source("butils")
# library(butils)
library(testthat)
set.seed(10)

context("#### R2 ####")

library(lava)
library(nlme)

m <- lvm(Y~X1+X2+X3)
regression(m) <- G~1
categorical(m, K=5, label = c("A","B","C","D","E")) <- ~G
d <- lava::sim(m, 1e2)
d <- d[order(d$G),]
d$Y <- d$Y + 0.5*as.numeric(as.factor(d$G))

## linear model
m.lm <- lm(Y~X1+X2+X3, data = d)
print(m.lm)
print(d)
res1 <- calcR2(m.lm)
res2 <- heplots::etasq(m.lm)

print(res1)

expect_equal(as.double(summary(m.lm)$r.squared), as.double(res1$R2.Buse))
expect_equal(as.double(summary(m.lm)$r.squared), as.double(res1$R2.McFadden))

expect_equal(as.double(na.omit(res2$`Partial eta^2`)), 
             as.double(res1$partialR2.Buse))
expect_equal(as.double(na.omit(res2$`Partial eta^2`)),
             as.double(res1$partialR2.McFadden))

