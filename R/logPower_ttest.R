### logPower_ttest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2020 (18:30) 
## Version: 
## Last-Updated: dec 17 2020 (21:32) 
##           By: Brice Ozenne
##     Update #: 219
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
##' @param method.meanvar [character] method used to identify the expected mean and variance on the log scale.
##' Either assumes that the outcome is log-normally distributed on the original scale (\code{method.meanvar = "lognorm"}),
##' or that it is normally distributed on the log-scale (\code{method.meanvar = "lognorm2"}).
##' @param method.cor [character] method used to identify the correlation on the log scale: one of \code{"uniroot"}, \code{"optim"}, or \code{"taylor"}.
##' @param sig.level [numeric 0-1] type 1 error
##' @param power [numeric 0-1] statistical power (i.e. complement to 1 of the type 2 error)
##' @param type [character] type of study:
##' difference in means between two independent groups (\code{type="two.sample"}),
##' difference in means between paired measurements (\code{type="paired"}),
##' equivalence in means between paired measurements (\code{type="equivalence"}),
##' @param n.large [integer, >0] sample size used to indentify the correlation coefficient or assess the error made when identifying the parameters. Should be large.
##' @param trace [logical] Should a progress bar be displayed when estimating the power for various combinaisons of parameters?
##' @param ncpus [integer, >0] Number of cores to be used, i.e., how many processes can be run simultaneously.
##'
##' @references
##' Shein-Chung Chow , Jun Shao & Hansheng Wang (2002). A note on sample size calculation for mean comparisons based on noncentral t-statistics, Journal of Biopharmaceutical Statistics, 12:4, 441-456, DOI: 10.1081/BIP-120016229


## * logPower_ttest - examples
##' @rdname logPower_ttest
##' @examples
##' if(require(MESS) && require(mvtnorm)){
##'
##' #### two sample comparison: 30% increase ####
##' X <- rlnorm(1e5, meanlog = 0, sdlog = 0.5)
##' Y <- rlnorm(1e5, meanlog = 0.275, sdlog = 0.5)
##'
##' ## proposed solution
##' logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
##'                gamma = 0.3, type = "two.sample")
##' 
##' ## no log-transform
##' beta <- mean(X)*0.3
##' sigma <- var(X)
##' power_t_test(delta = beta/sigma, power = 0.80)
##' 
##' ## using log-transform
##' beta <- mean(log(Y))-mean(log(X))
##' sigma <- sqrt(var(log(Y))/2+var(log(X))/2)
##' power_t_test(delta = beta/sigma, power = 0.80)
##'
##' ## using simulation
##' warper <- function(i, n){
##' X <- rlnorm(n, meanlog = 0, sdlog = 0.5)
##' Y <- rlnorm(n, meanlog = 0.275, sdlog = 0.5)
##' t.test(log(X),log(Y))$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 53))<=0.05)
##'
##' ## graphical display
##' df.power <- logPower_ttest(mu.original = mean(X),
##'                sigma2.original = var(X),
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
##' #### one sample comparison: 30% increase ####
##' rho <- 0.3
##' Sigma <- matrix(c(0.5^2,0.5^2*rho,0.5^2*rho,0.5^2),2,2)
##' XY <- exp(rmvnorm(1e5, mean = c(0,0.275), sigma = Sigma))
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
##' ## no log-transform
##' beta <- mean(X)*0.3
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
##' XY <- rmvnorm(n, mean = c(0,0.275), sigma = Sigma)
##' t.test(XY[,2]-XY[,1])$p.value
##' }
##' mean(unlist(lapply(1:1e3,warper, n = 38))<=0.05)
##' }
##' 
##' #### one sample equivalence: 10% difference ####
##' mu <- rep(0,2)
##' rho <- 0.3
##' Sigma <- matrix(c(0.5^2,0.5^2*rho,0.5^2*rho,0.5^2),2,2)
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
##' ls.res <- do.call(rbind,lapply(1:1e3,warper, n = 300))
##' mean(ls.res[,4]<=0.05)
##' }

## * logPower_ttest - code
##' @rdname logPower_ttest
##' @export
logPower_ttest <- function(mu.original, sigma2.original, gamma, 
                           rho = NULL, n = NULL, ratio.n = 1,
                           sig.level = 0.05, power = 0.8,
                           type = "two.sample", equivalence = FALSE, 
                           method.meanvar = "lognorm", method.cor = "uniroot",
                           n.large = 1e5, trace = TRUE, ncpus = NULL){

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
        attr(n.large,"diagnostic") <- FALSE

        warper <- function(iG){
              iOut <- logPower_ttest(mu.original = grid[iG,"mu.original"],
                                   sigma2.original = grid[iG,"sigma2.original"],
                                   gamma = grid[iG,"gamma"], 
                                   rho = if("rho" %in% name.grid){grid[iG,"rho"]}else{NULL},
                                   n = if("n" %in% name.grid){grid[iG,"n"]}else{NULL},
                                   ratio.n = grid[iG,"ratio.n"],
                                   method.meanvar = method.meanvar,
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
    type <- match.arg(type, c("two.sample","paired"))
    method.cor <- match.arg(method.cor, c("taylor","uniroot","optim"))
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

    ## ** identify mean and variance on the log-scale
    if(method.meanvar=="lognorm"){
        s0 <- log(1+sigma2.original/mu.original^2)
        a0 <- log(mu.original) - s0/2
        ## s1 <- log(1+sigma2.original/(mu.original*(1+gamma))^2)
        ## a1 <- log(mu.original*(1+gamma)) - s1/2
    }else if(method.meanvar=="lognorm2"){
        s0 <- uniroot(function(x){
            mu.original^2/sigma2.original - (1+x/2+x^2/8+x^3/48)^2/(x+(3/2)*x^2+(7/6)*x^3+(11/24)*x^4+(21/320)*x^5)},
            interval = c(1e-12,sigma2.original))$root
        a0 <- log(mu.original) - log(1+s0/2+s0^2/8+s0^3/48)
    }
    a1 <- a0 + log(1+gamma)

    ## ** identify correlation and pooled variance
    if(type %in% "paired"){
        
        if(method.cor=="taylor"){
            rho.log <- uniroot(function(x){
                rho - (x+1.5*x^2*s0+(1/12)*s0^2*(2*x^3+3*x))/(1+(3/2)*s0+(7/6)*s0^2+(11/24)*s0^3+(21/320)*s0^4)
            },interval = c(0,0.9999))$root
        }else if(method.cor=="uniroot"){
            require(mvtnorm)

            rho.log <- uniroot(f = function(x){
                Sigma.log <- matrix(c(s0,x*s0,x*s0,s0),2,2)
                return(cor(exp(mvtnorm::rmvnorm(n.large,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho)
            }, lower = 0, upper = 0.999)$root ## make sure that Sigma.log is positive definite
        }else if(method.cor=="optim"){
            require(mvtnorm)

            rho.log <- optim(fn = function(x){
                Sigma.log <- matrix(c(s0,x*s0,x*s0,s0),2,2)
                diff <- cor(exp(mvtnorm::rmvnorm(n.large,mean = c(a0,a1), sigma = Sigma.log)))[1,2]-rho
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
    if(!identical(attr(n.large,"diagnostic"),FALSE)){
        require(mvtnorm)

        Sigma.log <- matrix(c(s0,rho.log*s0,rho.log*s0,s0),2,2)
        Z <- exp(mvtnorm::rmvnorm(n.large,mean = c(a0,a1), sigma = Sigma.log))

        df.diagnostic <- rbind(requested = data.frame(mu=mu.original,sigma2=sigma2.original,gamma=gamma,rho=rho),
                               error = data.frame(mu=mu.original-mean(Z[,1]),sigma2=sigma2.original-var(Z[,1]),gamma=gamma-(mean(Z[,2])/mean(Z[,1])-1),rho=rho-cor(Z[,1],Z[,2]))
                               )
    }else{
        df.diagnostic <- NULL
    }
    
    ## ** statistical test
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
        tt$alternative <- "two.sided"
        
    }else{
        require(MESS)
        tt <- MESS::power_t_test(n = n,
                                 delta = a1-a0,
                                 sd = sqrt(s.pool),
                                 power = power,
                                 sig.level = sig.level,
                                 type = type,
                                 ratio = ratio.n)

        if(type == "paired"){
            tt$note <- "n is the number of *pairs"
        }
    }

    ## ** export
    attr(tt,"param.log.normal") <- c("a0"=a0,
                                     "a1"=a1,
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
