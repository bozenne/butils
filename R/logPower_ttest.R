### logPower_ttest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2020 (18:30) 
## Version: 
## Last-Updated: dec  5 2020 (15:58) 
##           By: Brice Ozenne
##     Update #: 153
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
##' @param mu.original [numeric] expected value of the outcome on the original scale for the control group or baseline.
##' @param sigma2.original [numeric] variance of the outcome on the original scale. Assumed to be the same over groups or over time.
##' @param gamma [numeric] relative change for the expected value in the treatment group or follow-up value.
##' @param rho [numeric 0-1] correlation coefficient between the paired values on the original scale. Only used when \code{type="paired"} or \code{type="equivalence"}.
##' @param n [integer >1] Sample size.
##' @param ratio.n [numeric >1] ratio between the sample size in larger and the smaller group. By default 1.
##' Note that the larger group is the reference group (with expectation \code{mu.original}) and the smaller group the active group (with expectation \code{gamma*mu.original}).
##' @param method [character] method used to identify the expected mean and variance on the log scale.
##' Either assumes that the outcome is log-normally distributed on the original scale (\code{method = "lognorm"}),
##' or that it is normally distributed on the log-scale (\code{method = "lognorm"}).
##' @param method.search [character] method used to identify the correlation on the log scale: either \code{"uniroot"} or \code{"optim"}.
##' @param sig.level [numeric 0-1] type 1 error
##' @param power [numeric 0-1] statistical power (i.e. complement to 1 of the type 2 error)
##' @param type [character] type of study:
##' difference in means between two independent groups (\code{type="two.sample"}),
##' difference in means between paired measurements (\code{type="paired"}),
##' equivalence in means between paired measurements (\code{type="equivalence"}),
##'
##' @references
##' Shein-Chung Chow , Jun Shao & Hansheng Wang (2002) A NOTE ONSAMPLE SIZE CALCULATION FOR MEAN COMPARISONS BASED ON NONCENTRAL t-STATISTICS, Journal of Biopharmaceutical Statistics, 12:4, 441-456, DOI: 10.1081/BIP-120016229


## * logPower_ttest - examples
##' @rdname logPower_ttest
##' @examples
##' if(require(MESS) && require(mvtnorm)){
##'
##' #### two sample comparison: 30% increase ####
##' X <- rlnorm(1e5, meanlog = 0, sdlog = 1)
##' Y <- rlnorm(1e5, meanlog = 0.41, sdlog = sqrt(0.69))
##'
##' ## proposed solution
##' logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.3, type = "two.sample")
##' 
##' ## no log-transform
##' beta <- mean(Y)-mean(X)
##' sigma <- sqrt(var(Y)/2+var(X)/2)
##' power_t_test(delta = beta/sigma, power = 0.80)
##' 
##' ## using log-transform
##' beta <- mean(log(Y))-mean(log(X))
##' sigma <- sqrt(var(log(Y))/2+var(log(X))/2)
##' power_t_test(delta = beta/sigma, power = 0.80)
##'
##' ## using simulation
##' warper <- function(i, n){
##' X <- rlnorm(n, meanlog = 0, sdlog = 1)
##' Y <- rlnorm(n, meanlog = 0.41, sdlog = sqrt(0.69))
##' t.test(log(X),log(Y))$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 80))<=0.05)
##'
##' ## graphical display
##' df.power <- logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = seq(0.1,0.4, length.out = 30),
##'                n = seq(50, 100, length.out = 30),
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
##' }
##' 
##' #### one sample comparison: 30% increase ####
##' rho <- 0.3
##' XY <- exp(rmvnorm(1e5, mean = c(0,0.41), sigma = rho+diag(c(1,0.69)-rho,2,2)))
##' X <- XY[,1]
##' ## colMeans(XY)-c(mean(X),mean(Y))
##' ## apply(XY,2,sd)-c(sd(X),sd(Y))
##'
##' ## proposed solution
##' logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.3, rho = cor(XY[,1],XY[,2]),
##'                type = "paired")
##'
##' 
##' ## no log-transform
##' beta <- mean(XY[,2]-XY[,1])
##' sigma <- sd(XY[,2]-XY[,1])
##' power_t_test(delta = beta/sigma, power = 0.80, type = "one.sample")
##' 
##' ## using log-transform
##' beta <- mean(log(XY[,2])-log(XY[,1]))
##' sigma <- sd(log(XY[,2])-log(XY[,1]))
##' power_t_test(delta = beta/sigma, power = 0.80, type = "one.sample")
##'
##' ## using simulation
##' warper <- function(i, n){
##' XY <- rmvnorm(n, mean = c(0,0.41), sigma = rho+diag(c(1,0.69)-rho,2,2))
##' t.test(XY[,2]-XY[,1])$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 52))<=0.05)
##' }
##' 
##' #### one sample equivalence: 10% difference ####
##' mu <- rep(0.1,2)
##' rho <- 0.5
##' Sigma <- 0.1*(rho+diag(c(1,1)-rho,2,2))
##' gamma <- 0.1
##' XY <- exp(rmvnorm(1e4, mean = mu, sigma = Sigma))
##' X <- XY[,1]
##'
##' ## proposed solution
##' logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = gamma, rho = cor(XY[,1],XY[,2]), n = 100, 
##'                type = "paired", equivalence = TRUE)
##' 
##' ## using log-transform
##' power_t_test(delta = log(1 + gamma)/sd(log(XY[,2]/XY[,1])), power = 1 - (1-0.80)/2,
##'              type = "one.sample", alternative = "one.sided")
##' ##pwr::pwr.t.test(d = log(1 + gamma)/sd(log(XY[,2]/XY[,1])), power = 1 - (1-0.80)/2,
##' ##             type = "one.sample", alternative = "greater")
##' power_t_test(delta = -log(1 - gamma)/sd(log(XY[,2]/XY[,1])), power = 1 - (1-0.80)/2,
##'              type = "one.sample", alternative = "one.sided")
##' ##pwr::pwr.t.test(d = log(1 - gamma)/sd(log(XY[,2]/XY[,1])), power = 1 - (1-0.80)/2,
##' ##             type = "one.sample", alternative = "less")
##' 
##' ## using simulation
##' warper <- function(i, n){
##' XY <- exp(rmvnorm(n, mean = mu, sigma = Sigma))
##' Z <- log(XY[,2]/XY[,1])
##' tt1 <- t.test(Z, mu = log(0.9), alternative = "greater")$p.value
##' tt2 <- t.test(Z, mu = log(1.1), alternative = "less")$p.value
##' return(c(colMeans(XY),mean(Z),max(tt1,tt2)))
##' }
##' ls.res <- do.call(rbind,lapply(1:1e3,warper, n = 89))
##' mean(ls.res[,4]<=0.05)
##' }

## * logPower_ttest - code
##' @rdname logPower_ttest
##' @export
logPower_ttest <- function(mu.original, sigma2.original, gamma, 
                           rho = NULL, n = NULL, ratio.n = 1,
                           sig.level = 0.05, power = 0.8,
                           type = "two.sample", equivalence = FALSE, 
                           method = "lognorm", method.search = "uniroot"){
    require(MESS)
    require(mvtnorm)

    ## ** deal with vector case
    ll <- list(mu.original = mu.original,
               sigma2.original = sigma2.original,
               gamma = gamma,
               rho = rho,
               n = n,
               ratio.n = ratio.n)
    grid <- do.call("expand.grid",ll[sapply(ll,length)>0])
    n.grid <- NROW(grid)
    if(n.grid>1){
        name.grid <- names(grid)
        
        ls.power <- lapply(1:n.grid, function(iG){
            iOut <- logPower_ttest(mu.original = grid[iG,"mu.original"],
                                   sigma2.original = grid[iG,"sigma2.original"],
                                   gamma = grid[iG,"gamma"], 
                                   rho = if("rho" %in% name.grid){grid[iG,"rho"]}else{NULL},
                                   n = if("n" %in% name.grid){grid[iG,"n"]}else{NULL},
                                   ratio.n = grid[iG,"ratio.n"],
                                   method = method,
                                   sig.level = sig.level,
                                   power = if("power" %in% name.grid){grid[iG,"power"]}else{NULL},
                                   type = type,
                                   method.search = method.search)
            return(cbind(grid[iG,setdiff(name.grid,c("n","power"))],
                         power = iOut$power,
                         n1 = iOut$n[1],
                         n2 = iOut$n[2]))
        })
        df.out <- do.call(rbind,ls.power)
        return(df.out)
       
    }
    
    ## ** check and normalize arguments
    method <- match.arg(method, c("lognorm","lognorm2"))
    type <- match.arg(type, c("two.sample","paired"))
    method.search <- match.arg(method.search, c("uniroot","optim"))
    if(type %in% c("paired") && ratio.n!=1){
        stop("Argument \'ratio.n\' must be 1 when argument \'type\' is \"paired\". \n")
    }
    if(is.null(n) && ratio.n!=1){
        stop("Argument \'ratio.n\' must be 1 when argument \'n\' not specified. \n")
    }
    if(type %in% c("paired","equivalence") && is.null(rho)){
        stop("Argument \'rho\' must be specified when argument \'type\' is \"paired\". \n")
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

    if(type %in% "paired"){
        if(method.search=="uniroot"){
            rho.log <- uniroot(f = function(x){
                Sigma.log <- matrix(c(s0,x*sqrt(s0)*sqrt(s1),x*sqrt(s0)*sqrt(s1),s1),2,2)
                return(cor(exp(mvtnorm::rmvnorm(1e5,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho)
            }, lower = 0, upper = min(sqrt(s1/s0),sqrt(s0/s1)))$root ## make sure that Sigma.log is positive definite
        }else if(method.search=="optim"){
            rho.log <- optim(fn = function(x){
                Sigma.log <- matrix(c(s0,x*sqrt(s0)*sqrt(s1),x*sqrt(s0)*sqrt(s1),s1),2,2)
                diff <- cor(exp(rmvnorm(1e5,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho
                return(diff^2)
            }, par = rho, lower = 0, upper = min(sqrt(s1/s0),sqrt(s0/s1)), method = "L-BFGS-B")$par ## make sure that Sigma.log is positive definite
        }
        s.pool <- s0+s1-2*sqrt(s0)*sqrt(s1)*rho.log
    }else{
        if(!is.null(n)){
            s.pool <- (s0*(n*ratio.n)+s1*n)/(n*(1+ratio.n))
        }else{
            s.pool <- s0/2+s1/2
        }
    }

    if(equivalence){

        power2 <- 1-(1-power)/2
        if(type == "two.sample"){
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
        }else if(type == "paired"){
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
    }else{
        tt <- MESS::power_t_test(n = n,
                                 delta = a1-a0,
                                 sd = sqrt(s.pool),
                                 power = power,
                                 sig.level = sig.level,
                                 type = type,
                                 ratio = ratio.n)
    }
    return(tt)    
}

######################################################################
### logPower_ttest.R ends here
