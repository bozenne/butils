### logPower_ttest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2020 (18:30) 
## Version: 
## Last-Updated: jan  6 2021 (17:29) 
##           By: Brice Ozenne
##     Update #: 293
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logPower_ttest
##' @title Power Calculation For t-tests on Log-Transformed Data
##' @description Power calculation for t-tests on log-transformed data.
##' @name logPower_ttest
##' 
##' @param type [character] type of study (parameter of interest):
##' \code{type="paired one.sample"} (mean change within the group), or
##' \code{type="two.sample"} (difference in mean between two independent groups), or
##' \code{type="paired two.sample"} (difference in mean change between two independent groups)
##' @param equivalence If TRUE, the alternative hypothesis is no change/difference. If FALSE, the null hypothesis is no change/difference.
##' @param gamma [numeric] When considering two groups: relative difference in expected values between the two groups on the original scale.
##' When considering a single group: expected value of the individual relative differences on the original scale.
##' @param sigma2 [numeric] When the argument \code{mu} is unspecified: coefficient of variation (i.e. standard error divided by expectation) of the outcome on the original scale.
##' When the argument \code{mu} is specified: variance of the outcome on the original scale.
##' @param mu [numeric] expected value of the outcome on the original scale for the control group or at baseline.
##' @param rho [numeric 0-1] correlation coefficient between the paired values on the original scale. 
##' @param n [integer >1] Sample size.
##' @param ratio.n [numeric >1] ratio between the sample size in larger and the smaller group. By default 1.
##' Note that the larger group is the reference group (with expectation \code{mu}) and the smaller group the active group (with expectation \code{gamma*mu}).
##' @param sig.level [numeric 0-1] type 1 error
##' @param power [numeric 0-1] statistical power (i.e. complement to 1 of the type 2 error)
##' @param method.meanvar [character] method used to identify the expected mean and variance on the log scale.
##' Either assumes that the outcome is log-normally distributed on the original scale (\code{method.meanvar = "lognorm"}),
##' or that it is normally distributed on the log-scale (\code{method.meanvar = "lognorm2"}).
##' @param method.cor [character] method used to identify the correlation on the log scale: one of \code{"uniroot"}, \code{"optim"}, or \code{"taylor"}.
##' @param n.large [integer, >0] sample size used to indentify the correlation coefficient or assess the error made when identifying the parameters. Should be large.
##' @param trace [logical] Should a progress bar be displayed when estimating the power for various combinaisons of parameters?
##' @param ncpus [integer, >0] Number of cores to be used, i.e., how many processes can be run simultaneously.
##'
##' @details The power calculation is performed under the assumption that on the log-scale, the variance of the outcome is the same between the two groups or over time.
##' 
##' @references
##' Shein-Chung Chow , Jun Shao & Hansheng Wang (2002). A note on sample size calculation for mean comparisons based on noncentral t-statistics, Journal of Biopharmaceutical Statistics, 12:4, 441-456, DOI: 10.1081/BIP-120016229
##' Gerald Van Belle and Donald C. Martin. Sample size as a function of coefficient of variation and ratio of means (1993). The American Statistician 47(3) 165-167.


## * logPower_ttest - examples
##' @examples
##' if(require(MESS) && require(mvtnorm)){
##'
##' rho <- 0.6
##' mu <- c(0,0.275)
##' Sigma <- matrix(c(0.5^2,0.5^2*rho,0.5^2*rho,0.5^2),2,2)
##'
##' #### 1- one sample with paired data ####
##' X <- exp(rmvnorm(1e5, mean = mu, sigma = Sigma))
##'
##' ## 1a) H0: no difference
##' # proposed solution
##' logPower_ttest(mu = mean(X[,1]),
##'                sigma2 = var(X[,1]),
##'                gamma = 0.3, rho = cor(X[,1],X[,2]),
##'                type = "paired one.sample")
##' 
##' logPower_ttest(sigma2 = sd(X[,1]) / mean(X[,1]),
##'                gamma = 0.3, rho = cor(X[,1],X[,2]),
##'                type = "paired one.sample")
##' 
##' # no log-transform
##' beta <- mean(X[,1])*0.3
##' sigma <- sd(X[,2]-X[,1])
##' power_t_test(delta = beta/sigma, power = 0.80, type = "one.sample")
##' 
##' # using log-transform
##' beta <- mean(log(X[,2])-log(X[,1]))
##' sigma <- sd(log(X[,2])-log(X[,1]))
##' power_t_test(delta = beta/sigma, power = 0.80, type = "one.sample")
##'
##' # using simulation
##' warper <- function(i, n){
##' X <- rmvnorm(n, mean = mu, sigma = Sigma)
##' t.test(X[,2]-X[,1])$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 39))<=0.05)
##' 
##' ## 1b) H0: >10% difference
##' # proposed solution
##' logPower_ttest(mu = mean(X[,1]),
##'                sigma2 = var(X[,1]),
##'                gamma = 0.1, rho = cor(X[,1],X[,2]),
##'                type = "paired one.sample", equivalence = TRUE)
##' 
##' ## using log-transform
##' power_t_test(delta = log(1 + 0.1)/sd(log(X[,2]/X[,1])), power = 1 - (1-0.80)/2,
##'              type = "one.sample", alternative = "one.sided")
##' power_t_test(delta = -log(1 - 0.1)/sd(log(X[,2]/X[,1])), power = 1 - (1-0.80)/2,
##'              type = "one.sample", alternative = "one.sided")
##' 
##' ## using simulation
##' warper <- function(i, n){
##' X <- exp(rmvnorm(n, mean = c(0,0), sigma = Sigma))
##' Z <- log(X[,2]/X[,1])
##' tt1 <- t.test(Z, mu = log(0.9), alternative = "greater")$p.value
##' tt2 <- t.test(Z, mu = log(1.1), alternative = "less")$p.value
##' return(c(colMeans(X),mean(Z),max(tt1,tt2)))
##' }
##' ls.res <- do.call(rbind,lapply(1:1e3,warper, n = 305))
##' mean(ls.res[,4]<=0.05)
##'
##' #### 2- two independent samples ####
##' X <- rlnorm(1e5, meanlog = mu[1], sdlog = 0.5)
##' Y <- rlnorm(1e5, meanlog = mu[2], sdlog = 0.5)
##'
##' ## 2a) H0: no difference
##' # proposed solution
##' logPower_ttest(mu = mean(X),
##'                sigma2 = var(X),
##'                gamma = 0.3, type = "two.sample")
##' 
##' logPower_ttest(sigma2 = sd(X)/mean(X),
##'                gamma = 0.3, type = "two.sample")
##' 
##' ## no log-transform
##' beta <- mean(X)*0.3
##' sigma <- sd(X)
##' power_t_test(delta = beta/sigma, power = 0.80, type = "two.sample")
##' 
##' ## using log-transform
##' beta <- mean(log(Y))-mean(log(X))
##' sigma <- sqrt(var(log(Y))/2+var(log(X))/2)
##' power_t_test(delta = beta/sigma, power = 0.80, type = "two.sample")
##'
##' ## using simulation
##' warper <- function(i, n){
##' X <- rlnorm(n, meanlog = mu[1], sdlog = 0.5)
##' Y <- rlnorm(n, meanlog = mu[2], sdlog = 0.5)
##' t.test(log(X),log(Y))$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 53))<=0.05)
##'
##' ## graphical display
##' df.power <- logPower_ttest(mu = mean(X),
##'                sigma2 = var(X),
##'                gamma = seq(0.1,0.4, length.out = 30),
##'                n = seq(20, 75, length.out = 30),
##'                type = "two.sample")
##'
##' if(require(ggplot2) && require(metR)){
##' gg <- ggplot(df.power, aes(x = n1, y = gamma))
##' gg <- gg + geom_tile(aes(fill = power))
##' gg <- gg + geom_contour(aes(z = power))
##' gg <- gg + scale_fill_gradientn(colours = terrain.colors(10), limits = c(0,1))
##' gg <- gg + geom_label_contour(aes(z = power), skip = 0)
##' gg <- gg +  scale_y_continuous(labels = scales::percent)
##' gg <- gg + xlab("sample size in each group") + ylab("mean difference between the groups")
##' gg
##' }
##' 
##' #### 3- two independent samples with paired data ####
##' X <- exp(rmvnorm(1e5, mean = rep(mu[1],2), sigma = Sigma))
##' Y <- exp(rmvnorm(1e5, mean = mu, sigma = Sigma))
##' 
##' ## 3a) H0: no difference
##' # proposed solution
##' logPower_ttest(mu = mean(X[,1]),
##'                sigma2 = var(X[,1]),
##'                gamma = 0.3, rho = cor(X[,1],X[,2]),
##'                type = "paired two.sample")
##' 
##' logPower_ttest(sigma2 = sd(X[,1])/mean(X[,1]),
##'                gamma = 0.3, rho = cor(X[,1],X[,2]),
##'                type = "paired two.sample")
##' 
##' ## no log-transform
##' beta <- mean(X[,1])*0.3
##' sigma <- sd(X[,1])
##' power_t_test(delta = beta/sigma, power = 0.80, type = "two.sample")
##'
##' ## using log-transform
##' beta <- mean(log(Y[,2])-log(Y[,1]))-mean(log(X[,2])-log(X[,1]))
##' sigma <- sqrt(var(log(Y[,2])-log(Y[,1]))/2+var(log(X[,2])-log(X[,1]))/2)
##' power_t_test(delta = beta/sigma, power = 0.80, type = "two.sample")
##'
##' ## using simulation
##' warper <- function(i, n){
##' X <- exp(rmvnorm(n, mean = rep(mu[1],2), sigma = Sigma))
##' Y <- exp(rmvnorm(n, mean = mu, sigma = Sigma))
##' t.test(log(X[,2])-log(X[,1]),log(Y[,2])-log(Y[,1]))$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 73))<=0.05)
##' }

## * logPower_ttest - code
##' @export
logPower_ttest <- function(type, equivalence = FALSE,
                           gamma, sigma2, mu = NULL, rho = NULL, n = NULL, ratio.n = 1,
                           sig.level = 0.05, power = 0.8,
                           method.meanvar = "lognorm", method.cor = "uniroot",
                           n.large = 1e5, trace = TRUE, ncpus = NULL){

    ## ** deal with vector case
    ll <- list(mu = mu,
               sigma2 = sigma2,
               gamma = gamma,
               rho = rho,
               n = n,
               ratio.n = ratio.n)
    grid <- do.call("expand.grid",ll[sapply(ll,length)>0])
    n.grid <- NROW(grid)
    if(n.grid>1){
        name.grid <- names(grid)
        attr(n.large,"diagnostic") <- FALSE
        warper <- function(iG){
            iOut <- logPower_ttest(sigma2 = sqrt(grid[iG,"sigma2"])/grid[iG,"mu"],
                                   gamma = grid[iG,"gamma"], 
                                   rho = if("rho" %in% name.grid){grid[iG,"rho"]}else{NULL},
                                   n = if("n" %in% name.grid){grid[iG,"n"]}else{NULL},
                                   ratio.n = grid[iG,"ratio.n"],
                                   method.meanvar = "lognorm",
                                   sig.level = sig.level,
                                   power = if("power" %in% name.grid){grid[iG,"power"]}else{NULL},
                                   type = type,
                                   method.cor = method.cor,
                                   n.large = n.large)
            return(cbind(grid[iG,setdiff(name.grid,c("n","power"))],
                         power = iOut$power,
                         n1 = iOut$n[1],
                         n2 = iOut$n[2]))
        }

        if(trace){
            require(pbapply)
            ls.power <- pbapply::pblapply(1:n.grid, warper, cl = ncpus)
        }else if(!is.null(ncpus)){
            require(parallel)
            ls.power <- parallel::mclapply(1:n.grid, warper, mc.cores = ncpus)
        }else{
            ls.power <- lapply(1:n.grid, warper)
        }
        df.out <- do.call(rbind,ls.power)
        return(df.out)
       
    }
    
    ## ** check and normalize arguments
    method.meanvar <- match.arg(method.meanvar, c("lognorm","lognorm2"))
    if(is.null(mu) && method.meanvar=="lognorm2"){
        stop("Argument \'method.meanvar\' must be \"lognorm\" when argument \'mu\' is NULL. \n")
    }
    type <- match.arg(type, c("paired one.sample", "two.sample", "paired two.sample"))
    paired <- grepl("paired",type, fixed = TRUE)
    type2 <- trimws(gsub("paired","",type, fixed = TRUE), which = "both")

    method.cor <- match.arg(method.cor, c("taylor","uniroot","optim"))
    if(type == "paired one.sample" && ratio.n!=1){
        stop("Argument \'ratio.n\' must be 1 when considering one sample paired observations. \n")
    }
    if(equivalence && ratio.n!=1){
        stop("Argument \'ratio.n\' must be 1 when argument \'equivalence\' is TRUE. \n")
    }
    if(paired && is.null(rho)){
        stop("Argument \'rho\' must be specified when considering paired observations. \n")
    }
    if(equivalence){
        if(is.null(power)){
            stop("Argument \'power\' must be specified when argument \'equivalence\' is TRUE. \n")
        }
        if(is.null(sig.level)){
            stop("Argument \'sig.level\' must be specified when argument \'equivalence\' is TRUE. \n")
        }
        if(!is.null(n)){
            n <- NULL
        }
    }

    ## ** identify mean and variance on the log-scale
    if(!is.null(mu)){ 
        if(method.meanvar=="lognorm"){
            s0 <- log(1+sigma2/mu^2)
            m0 <- log(mu) - s0/2
        }else if(method.meanvar=="lognorm2"){
            s0 <- uniroot(function(x){
                mu^2/sigma2 - (1+x/2+x^2/8+x^3/48)^2/(x+(3/2)*x^2+(7/6)*x^3+(11/24)*x^4+(21/320)*x^5)},
                interval = c(1e-12,sigma2))$root
            m0 <- log(mu) - log(1+s0/2+s0^2/8+s0^3/48)
        }
    }else{ 
        if(method.meanvar=="lognorm"){
            ## mu = exp(m0+s0/2)
            ## sigma2 = exp(2*m0+s0)(exp(s0-1))
            ## sigma2/mu^2 = exp(s0-1)
            s0 = log(1 + sigma2^2)
        }
        m0 <- 0 ## should not matter
    }
    
    if(equivalence){
        m1 <- m0
    }else{
        m1 <- m0 + log(1+gamma)
    }

    ## ** identify correlation and pooled variance
    if(paired){
        
        if(method.cor=="taylor"){
            rho.log <- uniroot(function(x){
                rho - (x+1.5*x^2*s0+(1/12)*s0^2*(2*x^3+3*x))/(1+(3/2)*s0+(7/6)*s0^2+(11/24)*s0^3+(21/320)*s0^4)
            },interval = c(0,0.9999))$root
        }else if(method.cor=="uniroot"){
            require(mvtnorm)

            rho.log <- uniroot(f = function(x){
                Sigma.log <- matrix(c(s0,x*s0,x*s0,s0),2,2)
                return(cor(exp(mvtnorm::rmvnorm(n.large,mean = c(m0,m1), sigma = Sigma.log)))[1,2]-rho)
            }, lower = 0, upper = 0.999)$root ## make sure that Sigma.log is positive definite
        }else if(method.cor=="optim"){
            require(mvtnorm)

            rho.log <- optim(fn = function(x){
                Sigma.log <- matrix(c(s0,x*s0,x*s0,s0),2,2)
                diff <- cor(exp(mvtnorm::rmvnorm(n.large,mean = c(m0,m1), sigma = Sigma.log)))[1,2]-rho
                return(diff^2)
            }, par = rho, lower = 0, upper = 0.999, method = "L-BFGS-B")$par ## make sure that Sigma.log is positive definite
        }
        s.pool <- 2*s0*(1-rho.log)
    }else{
        s.pool <- s0
        rho.log <- 0
        rho <- 0
    }

    ## ** diagnostic
    if(!identical(attr(n.large,"diagnostic"), FALSE)){
        require(mvtnorm)

        Sigma.log <- matrix(c(s0,rho.log*s0,rho.log*s0,s0),2,2)
        Z <- exp(mvtnorm::rmvnorm(n.large,mean = c(m0,m1), sigma = Sigma.log))

        if(!is.null(mu)){
            df.diagnostic <- rbind(requested = data.frame(mu=mu,sigma2=sigma2,gamma=gamma,rho=rho),
                                   error = data.frame(mu=mu-mean(Z[,1]),sigma2=sigma2-var(Z[,1]),gamma=gamma-(mean(Z[,2])/mean(Z[,1])-1),rho=rho-cor(Z[,1],Z[,2]))
                                   )
        }else{
            df.diagnostic <- rbind(requested = data.frame(sigma2=sigma2,gamma=gamma,rho=rho),
                                   error = data.frame(sigma2=sigma2-sd(Z[,1])/mean(Z[,1]),gamma=gamma-(mean(Z[,2])/mean(Z[,1])-1),rho=rho-cor(Z[,1],Z[,2]))
                                   )
        }
    }else{
        df.diagnostic <- NULL
    }
    
    ## ** statistical test
    if(equivalence){

        power2 <- 1-(1-power)/2
        if(type2 == "two.sample"){
            ## two sample equivalence test as in Shein-Chung Chow (2002)
            ## formula: n = 2 * (z_\alpha+z_(\beta/2))^2 / \delta^2
            tt <- list(pwr::pwr.t.test(d = log(1+gamma)/sqrt(s.pool),
                                       power = power2,
                                       sig.level = sig.level,
                                       type = "two.sample", alternative = "greater"),
                       pwr::pwr.t.test(d = log(1-gamma)/sqrt(s.pool),
                                       power = power2,
                                       sig.level = sig.level,
                                       type = "two.sample", alternative = "less"))
        }else if(type2 == "one.sample"){
            ## one sample equivalence test
            ## formula: n = (z_\alpha+z_(\beta/2))^2 / \delta^2
            tt <- list(pwr::pwr.t.test(d = log(1+gamma)/sqrt(s.pool),
                                       power = power2,
                                       sig.level = sig.level,
                                       type = "one.sample", alternative = "greater"),
                       pwr::pwr.t.test(d = log(1-gamma)/sqrt(s.pool),
                                       power = power2,
                                       sig.level = sig.level,
                                       type = "one.sample", alternative = "less"))
        }
        
        tt <- tt[[which.max(sapply(tt,"[[","n"))]]
        tt$power <- power
        tt$d <- log(1+gamma)/sqrt(s.pool)
        tt$method <- "Equivalence t test power calculation"
        tt$alternative <- "two.sided"
        
    }else{
        require(MESS)
        tt <- MESS::power_t_test(n = n,
                                 delta = log(1+gamma), ## same for all parametrisations
                                 sd = sqrt(s.pool),
                                 power = power,
                                 sig.level = sig.level,
                                 type = type2,
                                 ratio = ratio.n)

        if(type2 == "two.sample"){
            tt$note <- "n is the number of observation per group"
        }else if(paired){
            tt$note <- "n is the number of pairs"
        }
    }

    ## ** export
    attr(tt,"param.log.normal") <- c("m0"=m0,
                                     "m1"=m1,
                                     "s0"=s0,
                                     "s1"=s0,
                                     "rho.log"=rho.log)
    attr(tt,"diagnostic") <- df.diagnostic
    class(tt) <- append("logPower.htest",class(tt))
    return(tt)    
}

## * print.logPower.htest
print.logPower.htest <- function(x,...){

    param.log.normal <- attr(x,"param.log.normal")
    diagnostic <- attr(x,"diagnostic")

    class(x) <- setdiff(class(x), "logPower.htest")

    print(x)
    if(!is.null(diagnostic)){
        cat("      quality of the approximation:\n")
        print(diagnostic)
    }
    cat("\n")
    cat("      parametrisation of the log-normal distribution:\n")
    print(param.log.normal)
}

######################################################################
### logPower_ttest.R ends here
