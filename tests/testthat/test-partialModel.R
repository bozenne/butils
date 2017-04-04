library(testthat)

context("#### partialModel ####")

library(lava)
library(nlme)

m <- lvm(Y~X1+X2+X3+X4+X5)
categorical(m, K=3, labels = 1:3) <- ~X1
categorical(m, K=2, labels = letters[1:2]) <- ~X2
d <- sim(m, 1e2)


# {{{ lm
test_that("lm", {
  lm.full <- lm(Y~X1*X2+X3+X4+X5, data = d)
  
  res <- partialModel(lm.full, var1 = "X4", var2 = "X5")
  expect_equal(names(coef(res)), c("(Intercept)","X4","X5"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0)
  expect_equal(coef(res)[c("X4","X5")],coef(lm.full)[c("X4","X5")])
  
  res <- partialModel(lm.full, var1 = "X1", var2 = "X5")
  expect_equal(names(coef(res)), c("(Intercept)","X12","X13","X5"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0)
  expect_equal(coef(res)[c("X12","X13","X5")],coef(lm.full)[c("X12","X13","X5")])
  
  res <- partialModel(lm.full, var1 = "X2", var2 = "X3")
  expect_equal(names(coef(res)), c("(Intercept)","X2b","X3"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0)
  expect_equal(coef(res)[c("X2b","X3")],coef(lm.full)[c("X2b","X3")])
  
  res <- partialModel(lm.full, var1 = "X1", var2 = "X2")
  expect_equal(names(coef(res)), c("(Intercept)","X12","X13","X2b","X12:X2b","X13:X2b"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0)
  expect_equal(coef(res)[c("X12","X13","X2b","X12:X2b","X13:X2b")],coef(lm.full)[c("X12","X13","X2b","X12:X2b","X13:X2b")])
})
# }}}

# {{{ gls
test_that("gls", {
  gls.full <- gls(Y~X1*X2+X3+X4+X5, data = d,
                  weights = varIdent(form = ~ 1|X2))
  
  res <- partialModel(gls.full, var1 = "X4", var2 = "X5")
  expect_equal(names(coef(res)), c("(Intercept)","X4","X5"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(coef(res)[c("X4","X5")],coef(gls.full)[c("X4","X5")], tolerance = 1e-3)
  expect_equal(vcov(res)[c("X4","X5"),c("X4","X5")],  vcov(gls.full)[c("X4","X5"),c("X4","X5")],
               tolerance = 1e-2)
  
  res <- partialModel(gls.full, var1 = "X1", var2 = "X5")
  expect_equal(names(coef(res)), c("(Intercept)","X12","X13","X5"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(coef(res)[c("X12","X13","X5")],coef(gls.full)[c("X12","X13","X5")], tolerance = 1e-3)
  
  res <- partialModel(gls.full, var1 = "X2", var2 = "X3")
  expect_equal(names(coef(res)), c("(Intercept)","X2b","X3"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(coef(res)[c("X2b","X3")],coef(gls.full)[c("X2b","X3")], tolerance = 1e-3)
  
  res <- partialModel(gls.full, var1 = "X1", var2 = "X2")
  expect_equal(names(coef(res)), c("(Intercept)","X12","X13","X2b","X12:X2b","X13:X2b"))
  expect_equal(as.double(coef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(coef(res)[c("X12","X13","X2b","X12:X2b","X13:X2b")],coef(gls.full)[c("X12","X13","X2b","X12:X2b","X13:X2b")], tolerance = 1e-3)
})
# }}}

# {{{ lme
test_that("lme", {
  lme.full <- lme(Y~X1*X3+X4+X5, data = d, random = ~ 1|X2)
  
  res <- partialModel(lme.full, var1 = "X4", var2 = "X5")
  expect_equal(names(fixef(res)), c("(Intercept)","X4","X5"))
  expect_equal(as.double(fixef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(fixef(res)[c("X4","X5")],fixef(lme.full)[c("X4","X5")], tolerance = 1e-3)
  expect_equal(vcov(res)[c("X4","X5"),c("X4","X5")],  vcov(lme.full)[c("X4","X5"),c("X4","X5")],
               tolerance = 1e-2)
  
  res <- partialModel(lme.full, var1 = "X1", var2 = "X5")
  expect_equal(names(fixef(res)), c("(Intercept)","X12","X13","X5"))
  expect_equal(as.double(fixef(res)["(Intercept)"]),0, tolerance = 1e-3)
  expect_equal(fixef(res)[c("X12","X13","X5")],fixef(lme.full)[c("X12","X13","X5")], tolerance = 1e-3)
})
# }}}
