### spca.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 30 2021 (14:49) 
## Version: 
## Last-Updated: jul 30 2021 (18:14) 
##           By: Brice Ozenne
##     Update #: 39
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @references Yunjin Choi, Jonathan Taylor, and Robert Tibshirani. Selecting the number of principal components estimation of the true rank of a noisy matrix. The Annals of Statistics (2017) 45(6):2590-2617.
##' Joni Virta and Klaus Nordhausen. Estimating the number of signals using principal component analysis. Stat (2019).
##' 
##' @examples
##' set.seed(10)
##' n.obs <- 100
##' Z1 <- rnorm(n.obs, 0, 1)
##' Z2 <- rnorm(n.obs, 0, 1)
##' Z3 <- rnorm(n.obs, 0, 1)
##' X <- cbind(Z1 + rnorm(n.obs,0,0.25), Z2 + rnorm(n.obs,0,0.25), Z3 + rnorm(n.obs,0,0.25), Z1 + Z2 + rnorm(n.obs,0,0.25), Z1 - Z2 + rnorm(n.obs,0,0.25), Z1 - Z3 + rnorm(n.obs,0,0.25))
##'
##' princomp(X)
##' spca(X, method = "virta")
##'
##'
##' library(ICtest)
##' PCAasymp(X,k=0)
##' PCAasymp(X,k=1)
##' PCAasymp(X,k=2)
##' PCAasymp(X,k=3)
##' PCAasymp(X,k=4)
##' 
##' eigen.Sigma <- eigen(var(X))
##'
##' n.obs * NCOL(X) * var(eigen.Sigma$values)/(2*mean(eigen.Sigma$values)^4*1)
##' 
##' library(mvtnorm)
##' n.obs <- 1e3
##' sigma <- 10
##' p <- 10
##' beta <- 1
##' 
##' vec.stat <- sapply(1:1e3, function(iSim){
##' iSigma <- var(rmvnorm(n.obs, sigma = diag(sigma^2,p,p)))
##' n.obs * p * var(eigen(iSigma)$values)/(2*sigma^4*beta)
##' })
##'
##' library(ggplot2)
##' ggplot() + geom_histogram(data = data.frame(x=vec.stat), aes(x=x,y=..density..)) + geom_density(data = data.frame(x=rchisq(1e5, df = 0.5*(p-1)*(p+2))), aes(x=x), color = "red")
##' 

## * spca (code)
spca <- function(object, method, method.sigma2 = "mean"){

    ## ** check user input
    if(!is.matrix(object)){
        stop("Argument \'object\' must be a matrix. \n")
    }

    method <- match.arg(method, c("choi","virta"))
    if(!is.numeric(method.sigma2)){
        method.sigma2 <- match.arg(method.sigma2, c("mean"))
    }else{
        sigma2 <- method.sigma2
        method.sigma2 <- "manual"
    }

    ## ** prepare
    n.col <- NCOL(object)
    n.obs <- NROW(object)

    if(method == "choi"){
        objectN <- scale(object)
    }else{
        objectN <- object
    }
    Xc <- sweep(object, MARGIN = 2, FUN = "-", STATS = colMeans(X))
    Sigma <- crossprod(Xc)/n.obs
    beta <- mean(rowSums(Xc %*% solve(Sigma) * Xc)^2)/(n.col*(n.col+2))
    Sigma.eigen <- eigen(Sigma , symmetric = TRUE)
    lambda <- Sigma.eigen$values
    R <- Sigma.eigen$vector
    ## R %*% diag(lambda) %*% t(R) - Sigma

    out <- data.frame(matrix(NA,nrow=n.col, ncol=2, dimnames = list(NULL,c("statistic","p.value"))))
    
    ## ** spca
    for(iRank in 1:n.col){ ## iRank <- 1

        if(method.sigma2=="manual"){
            iSigma2.noise <- mean(lambda)
        }else if(method.sigma2=="mean"){
            iSigma2.noise <- mean(lambda[iRank:n.col])
        }
        iSigma2.lambda <-  sum((lambda[iRank:n.col]-mean(lambda[iRank:n.col]))^2)/(n.col-iRank+1)
        
    if(method == "choi"){
        lambda2 <- lambda^2
        bound <- c(5*sigma,lambda)
        calcIntegrand <- function(z,k){ ## z <- 4
            exp(-z^2/(2*sigma^2) + (n.obs-n.col)*log(z) + sum(abs(z^2 - lambda2[-k])))
        }
    
        calcRatio <- function(k){ ## k <- 1
            I1 <- integrate(f = function(z){calcIntegrand(z=z,k=k)}, lower = bound[k+1], upper = bound[k])
            I2 <- integrate(f = function(z){calcIntegrand(z=z,k=k)}, lower = bound[k+2], upper = bound[k])
            return(I1/I2)
        }
    }else if(method == "virta"){
        out$statistic[iRank] <- n.obs * (n.col-iRank+1) * iSigma2.lambda/(2*iSigma2.noise^2*beta)
        out$p.value[iRank] <- 1 - pchisq(out$statistic[iRank], df = 0.5*(n.col-iRank)*(n.col-iRank+3))
    }
    }

    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### spca.R ends here
