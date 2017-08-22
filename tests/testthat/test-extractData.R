context("#### Extract data from models ####")

set.seed(10)
n <- 100

#### linear regression ####
Y1 <- rnorm(n, mean = 0)
Y2 <- rnorm(n, mean = 0.3)
df <- rbind(data.frame(Y=Y1,G=1,Id = 1:5),
           data.frame(Y=Y2,G=2,Id = 1:5)
           )
m.lm <- lm(Y ~ G, data = df)
test_that("extractData (lm)", {
  expect_named(extractData(m.lm), expected = c("Y","G"))
})

library(nlme)
test_that("extractData (gls/lme/lmer)", {
  m.gls <- gls(Y ~ G, weights = varIdent(form = ~ 1|Id), data = df)
  expect_named(extractData(m.gls), expected = c("Y","G","Id"))
  m.lme <- lme(Y ~ G, random = ~ 1|Id, data = df)
  expect_named(extractData(m.lme), expected = c("Y","G","Id"))
  m.lmer <- lmer(Y ~ G + (1|Id), data = df)
  expect_named(extractData(m.lmer), expected = c("Y","G","Id"))
})


library(lava)
test_that("extractData (lvm)", {
  e <- estimate(lvm(Y ~ G), data = df)
  expect_named(extractData(e), expected = c("Y","G"))
})


#### survival ####
library(riskRegression)
library(survival)
dt <- sampleData(n, outcome = "survival")
test_that("extractData (survival)", {
  # no strata
  m.cox <- coxph(Surv(time, event) ~ X1 + X2, data = dt, x = TRUE, y = TRUE)
  expect_named(extractData(m.cox), expected = c("start","stop","status","X1","X2"))
  # strata
  m.cox <- coxph(Surv(time, event) ~ strata(X1) + X2, data = dt, x = TRUE, y = TRUE)
  expect_named(extractData(m.cox), expected = c("start","stop","status","X2","strata","X1"))
})
