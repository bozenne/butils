### logPower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2020 (18:30) 
## Version: 
## Last-Updated: dec  3 2020 (21:56) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * log.pwr.t.test
##' @title Power Calculation For t-tests on Log-Transformed Data
##' @description Power calculation for t-tests on log-transformed data.
##' 
##' @param mu.original [numeric] expected value of the outcome on the original scale for the control group or baseline.
##' @param sigma2.original [numeric] variance of the outcome on the original scale. Assumed to be the same over groups or over time.
##' @param gamma [numeric] relative change for the expected value in the treatment group or follow-up value.
##' @param rho [numeric 0-1] correlation coefficient between the paired values on the original scale. Only used when \code{type="paired"} or \code{type="equivalence"}.
##' @param method [character] method used to identify the expected mean and variance on the log scale.
##' Either assumes that the outcome is log-normally distributed on the original scale (\code{method = "lognorm"}),
##' or that it is normally distributed on the log-scale (\code{method = "lognorm"}).
##' @param sig.level [numeric 0-1] type 1 error
##' @param power [numeric 0-1] statistical power (i.e. complement to 1 of the type 2 error)
##' @param type [character] type of study:
##' difference in means between two independent groups (\code{type="two.sample"}),
##' difference in means between paired measurements (\code{type="paired"}),
##' equivalence in means between paired measurements (\code{type="equivalence"}),
##'
##' @examples
##' if(require(pwr) && require(mvtnorm)){
##'
##' #### two sample comparison: 30% increase ####
##' X <- rlnorm(1e5, meanlog = 0, sdlog = 1)
##' Y <- rlnorm(1e5, meanlog = 0.41, sdlog = sqrt(0.69))
##'
##' ## proposed solution
##' log.pwr.t.test(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.3, type = "two.sample")
##' 
##' ## no log-transform
##' beta <- mean(Y)-mean(X)
##' sigma <- sqrt(var(Y)/2+var(X)/2)
##' pwr.t.test(d = beta/sigma, power = 0.80)
##' 
##' ## using log-transform
##' beta <- mean(log(Y))-mean(log(X))
##' sigma <- sqrt(var(log(Y))/2+var(log(X))/2)
##' pwr.t.test(d = beta/sigma, power = 0.80)
##'
##' ## using simulation
##' warper <- function(i, n){
##' X <- rlnorm(n, meanlog = 0, sdlog = 1)
##' Y <- rlnorm(n, meanlog = 0.41, sdlog = sqrt(0.69))
##' t.test(log(X),log(Y))$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 80))<=0.05)
##' 
##' #### one sample comparison: 30% increase ####
##' rho <- 0.3
##' XY <- exp(rmvnorm(1e5, mean = c(0,0.41), sigma = rho+diag(c(1,0.69)-rho,2,2)))
##' ## colMeans(XY)-c(mean(X),mean(Y))
##' ## apply(XY,2,sd)-c(sd(X),sd(Y))
##'
##' ## proposed solution
##' log.pwr.t.test(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.3, rho = cor(XY[,1],XY[,2]),
##'                type = "paired")
##' 
##' ## no log-transform
##' beta <- mean(XY[,2]-XY[,1])
##' sigma <- sd(XY[,2]-XY[,1])
##' pwr.t.test(d = beta/sigma, power = 0.80, type = "one.sample")
##' 
##' ## using log-transform
##' beta <- mean(log(XY[,2])-log(XY[,1]))
##' sigma <- sd(log(XY[,2])-log(XY[,1]))
##' pwr.t.test(d = beta/sigma, power = 0.80, type = "one.sample")
##'
##' ## using simulation
##' warper <- function(i, n){
##' XY <- rmvnorm(n, mean = c(0,0.41), sigma = rho+diag(c(1,0.69)-rho,2,2))
##' t.test(XY[,2]-XY[,1])$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 52))<=0.05)
##' }
##' 
##' #### one sample equivalence: 15% difference ####
##' rho <- 0.3
##' XY <- exp(rmvnorm(1e5, mean = c(0,0), sigma = rho+diag(c(1,1)-rho,2,2)))
##'
##' ## proposed solution
##' log.pwr.t.test(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.15, rho = cor(XY[,1],XY[,2]),
##'                type = "equivalence")
##' 
##' ## using log-transform
##' beta <- mean(log(XY[,2])-log(XY[,1]))
##' sigma <- sd(log(XY[,2])-log(XY[,1]))
##' pwr.t.test(d = log(1.15)/(1+sigma^2/2), power = 1-(1-0.80)/2, type = "paired")
##'
##' ## using simulation
##' warper <- function(i, n){
##' XY <- rmvnorm(n, mean = c(0,0), sigma = rho+diag(c(1,1)-rho,2,2))
##' Z <- XY[,2]-XY[,1]
##' tt1 <- t.test(Z, mu = -log(1.15)/(1+var(Z)/2), alternative = "greater")$p.value
##' tt2 <- t.test(Z, mu = log(1.15)/(1+var(Z)/2), alternative = "less")$p.value
##' return(c(colMeans(exp(XY)),mean(exp(XY[,1]))/mean(exp(XY[,2])),max(tt1,tt2)))
##' }
##' ls.res <- do.call(rbind,lapply(1:1e3,warper, n = 1547))
##' mean(abs(ls.res[,3])>1.15)
##' hist(ls.res[,3])
##' mean(ls.res[,4]<=0.05)
##' }
##' @export
log.pwr.t.test <- function(mu.original, sigma2.original, gamma, rho,
                           method = "lognorm",
                           sig.level = 0.05,
                           power = 0.8,
                           type = c("two.sample")){
    require(pwr)
    require(mvtnorm)
    method <- match.arg(method, c("lognorm","lognorm2"))
    type <- match.arg(type, c("two.sample","paired","equivalence"))

    if(method=="lognorm"){
        s0 <- log(1+sigma2.original/mu.original^2)
        a0 <- log(mu.original) - s0/2
        s1 <- log(1+sigma2.original/(mu.original*(1+gamma))^2)
        a1 <- log(mu.original*(1+gamma)) - s1/2
    }else if(method=="lognorm2"){
        s0 <- uniroot(function(x){mu.original^2/sigma2.original - (1+x/2)^2/(x+x^2*7/4)},interval = c(0,1))$root
        a0 <- log(sigma2.original/(s0+s0^2*7/4))/2
        s1 <- uniroot(function(x){mu.original^2*(gamma+1)^2/sigma2.original - (1+x/2)^2/(x+x^2*7/4)},interval = c(0,1))$root
        a1 <- log(sigma2.original/(s1+s1^2*7/4))/2
    }

    if(type %in% c("paired","equivalence")){
        rho.log1 <- uniroot(f = function(rho.log){
            Sigma.log <- matrix(c(s0,rho.log*sqrt(s0)*sqrt(s1),rho.log*sqrt(s0)*sqrt(s1),s1),2,2)
            return(cor(exp(rmvnorm(1e5,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho)
        }, lower = 0, upper = min(sqrt(s1/s0),sqrt(s0/s1))) ## make sure that Sigma.log is positive definite

        ## rho.log2 <- optim(fn = function(rho.log){
        ##     Sigma.log <- matrix(c(s0,rho.log*sqrt(s0)*sqrt(s1),rho.log*sqrt(s0)*sqrt(s1),s1),2,2)
        ##     diff <- cor(exp(rmvnorm(1e5,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho
        ##     return(diff^2)
        ## }, par = rho.log1$root, lower = 0, upper = min(sqrt(s1/s0),sqrt(s0/s1)), method = "L-BFGS-B") ## make sure that Sigma.log is positive definite
        rho.log <- rho.log1$root
        s.pool <- s0+s1-2*sqrt(s0)*sqrt(s1)*rho.log
    }else{
        s.pool <- s0/2+s1/2
    }

    if(type=="equivalence"){
        warning("DO NOT USE: INCCORRECT RESULT!!!!!")
        return(pwr::pwr.t.test(d = log(1+gamma)/(1+s.pool/2), power = 1-(1-power)/2, sig.level = sig.level, type = "paired"))
    }else{
        return(pwr::pwr.t.test(d = (a1-a0)/sqrt(s.pool), power = power, sig.level = sig.level, type = "paired"))
    }
    
}

######################################################################
### logPower.R ends here
