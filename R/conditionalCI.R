### conditionalCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 20 2018 (14:43) 
## Version: 
## Last-Updated: nov 29 2018 (13:25) 
##           By: Brice Ozenne
##     Update #: 146
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * conditionalCI [documentation]
#' @title Conditional confidence interval 
#' @description Conditional confidence intervals for normally distributed random variables.
#' @name conditionalCI
#' 
#' @param theta [numeric vector] the observed value(s) of the statistic.
#' @param sigma [numeric vector] the standard error of the statistic(s).
#' @param threshold [numeric, >0] the threshold.
#' @param conf.level [numeric, 0-1] Conditional confidence level of the interval.
#' @param method [character] the method used to compute the conditional confidence interval.
#' So far only \code{"shortest"} is supported.
#' @param trace [interger, 0-2] should a progress bar be displayed?
#' @param ... additional arguments to be passed to \code{lavaSearch2:::.calcShortestCI}
#' to specify the optimization method.
#' 
#' @details Compute the confidence interval of theta
#' conditional on theta being greater (in absolute value) than the threshold.
#' theta is assumed to be normally distributed with known variance sigma.
#' The thresholding is applied after normalization (i.e. dividing theta by sigma).
#' This corresponding to conditioning on a minimal significance level.
#'
#' @examples
#' ci <- conditionalCI(theta = seq(3,7, by = 0.1), threshold = 2.9999)
#' print(ci)
#'
#' autoplot(ci)
#' 
#' ci1 <- conditionalCI(theta = 4, sigma = 2, threshold = 2.9999, method = "shortest")
#' ci2 <- conditionalCI(theta = 4/2, sigma = 1, threshold = 2.9999/2, method = "shortest")
#' confint(ci1) - 2*confint(ci2)
#' ci3 <- conditionalCI(theta = 4, sigma = 2, threshold = 2.9999, method = "shortest2")
#' confint(ci1) - confint(ci3)


## * conditionalCI [code]
#' @rdname conditionalCI
#' @export
conditionalCI <- function(theta, threshold,
                          sigma = 1, conf.level = 0.95, method = "shortest", 
                          trace = length(theta)>1, ...){

    ## ** normalize arguments
    n.theta <- length(theta)
    alpha <- 1 - conf.level
    method <- match.arg(method, choices = c("shortest","shortest2"))

    if(n.theta!=length(sigma)){
        if(length(sigma)==1){
            sigma <- rep(sigma, n.theta)
        }else{
            stop("The length of argument \'theta\' does not match the length of argument \'sigma\' \n",
                 "length(theta)=",n.theta,"\n",
                 "length(sigma)=",length(sigma),"\n")
        }            
    }
    if(length(threshold)!=1){
        stop("The length of argument \'threshold\' must be 1")
    }
    threshold <- rep(threshold, n.theta)

    if(any(threshold<0)){
        stop("Argument \'threshold\' must contain only positive values \n")
    }
    if(any(sigma<0)){
        stop("Argument \'sigma\' must contain only positive values \n")
    }
    ## ** reorder theta according to its relevance
    theta.save <- theta
    sigma.save <- sigma
    threshold.save <- threshold
    
    if(method == "shortest"){
        tempo.selected <- ((abs(theta)/sigma) > threshold/sigma)
        new.order <- order(  (abs(theta)/sigma) * sign(2*tempo.selected  - 1) )
    }else if(method == "shortest2"){
        new.order <- order(abs(theta))
    }
    
    theta <- theta[new.order]
    sigma <- sigma[new.order]
    threshold <- threshold[new.order]
    
    ## ** loop over theta
    M.AR <- matrix(NA, nrow = n.theta, ncol = 5,
                   dimnames = list(NULL, c("value","lower.left","lower.right","upper.left","upper.right")))
    M.CI <- matrix(NA, nrow = n.theta, ncol = 3,
                   dimnames = list(NULL, c("value","lower","upper")))
    df.optim <- NULL
    if(trace==1){pb <- utils::txtProgressBar(min = 0, max = n.theta, style = 3)}

    for(iSeq in 1:n.theta){
        if(method == "shortest"){
            iTheta <- as.double(theta[iSeq]/sigma[iSeq])
            iSigma <- 1
            iThreshold <- as.double(threshold[iSeq]/sigma[iSeq])
        }else if(method == "shortest2"){
            iTheta <- as.double(theta[iSeq])
            iSigma <- as.double(sigma[iSeq])
            iThreshold <- as.double(threshold[iSeq])
        }
        if(trace==1){utils::setTxtProgressBar(pb, iSeq)}        
        if(trace>1){cat(iTheta," ")}

        ## *** computation AR/CI
        if(abs(iTheta) < iThreshold){
            
            iRes <- list(AR = c(lower.left = NA, lower.right = NA, upper.left = NA, upper.right = NA),
                         CI = c(lower = NA, upper = NA),
                         optim = NULL)
              
        }else if(iThreshold == 0){
            
            iRes <- list(AR = c(lower.left = qnorm(alpha/2, mean = 0, sd = iSigma), lower.right = 0,
                                upper.left = 0, upper.right = qnorm(1 - alpha/2, mean = 0, sd = iSigma)),
                         CI = c(lower = iTheta + qnorm(alpha/2, mean = 0, sd = iSigma),
                                upper = iTheta + qnorm(1 - alpha/2, mean = 0, sd = iSigma)),
                         optim = NULL)

        }else{
            vec.init <- c(theta1 = abs(iTheta),
                          theta2 = abs(iTheta),
                          lower = abs(iTheta) + pnorm(alpha/2, mean = 0, sd = iSigma)*1,
                          upper = abs(iTheta) + pnorm(1 - alpha/2, mean = 0, sd = iSigma)*1)
            if(NROW(df.optim)>1){
                ## vec.init["lower"] <- ls.CI[[iSeq-1]]["lower"]
                vec.init["theta1"] <- df.optim[df.optim$param == "theta1" & df.optim$theta == theta[iSeq-1],"solution"]
                ## vec.init <- ls.CI[[iSeq-1]][c("theta1","theta2","lower","upper")]
            }
            iRes <- .calcShortestCI(value = abs(iTheta), sigma = iSigma, value.sign = sign(iTheta), threshold = iThreshold,
                                    alpha = alpha, vec.init = vec.init, ...)

            if(method == "shortest"){
                M.CI[iSeq,] <- c(theta[iSeq], iRes$CI * sigma[iSeq])
                M.AR[iSeq,] <- c(theta[iSeq], iRes$AR * sigma[iSeq])
            }else if(method == "shortest2"){
                M.CI[iSeq,] <- c(theta[iSeq], iRes$CI)
                M.AR[iSeq,] <- c(theta[iSeq], iRes$AR)
            }
        }

        ## store optim results
        if(!is.null(iRes$optim)){
            df.optim <- rbind(df.optim,
                              cbind(theta = theta[iSeq], iRes$optim)
                              )
        }
        
    
    }
    
    if(trace==1){close(pb)}
    if(trace>1){cat("\n")}
    
    ## ** export
    restaure.order <- order(new.order)
    out <- list(CI = M.CI[restaure.order,,drop=FALSE],
                AR = M.AR[restaure.order,,drop=FALSE],
                optim = df.optim,
                conf.level = conf.level,
                method = method,
                sigma = sigma[restaure.order],
                threshold = threshold[restaure.order])
    class(out) <- "conditionalCI"
    return(out)
}

## * print.conditionalCI
#' @export
print.conditionalCI <- function(x, ...){
    return(print(stats::confint(x)))
}

## * confint.conditionalCI
#' @export
confint.conditionalCI <- function(object, parm, level = 0.95, ...){
    return(as.data.frame(object$CI))
}

## * autoplot.conditionalCI
#' @title Display Conditional Confidence Intervals
#' @description Display conditional confidence intervals.
#'
#' @param object output of the function \code{conditionalCI}.
#' @param value [logical] should the identity line be displayed (in blue)?
#' @param unconditional.CI [logical] should the unconditional confidence intervals be displayed (in red)?
#' @param plot [logical] should the graph be displayed in a graphical window?
#' @param ... not used, for compatibility with the generic method.
#' 
#' @export
autoplot.conditionalCI <- function(object, value = TRUE, unconditional.CI = TRUE, plot = TRUE, ...){

    ## extract data
    CI <- stats::confint(object)

    ## create plot
    gg <- ggplot2::ggplot(CI, aes_string(x = "value"))
    gg <- gg + ggplot2::geom_ribbon(aes_string(ymin = "lower", ymax = "upper"))
    if(value){
        gg <- gg + ggplot2::geom_abline(intercept = 0, slope = 1, col = "blue")
    }
    if(unconditional.CI){
        gg <- gg + ggplot2::geom_abline(intercept = qnorm(0.975), slope = 1, col = "red")
        gg <- gg + ggplot2::geom_abline(intercept = qnorm(0.025), slope = 1, col = "red")
    }

    ## display
    if(plot){
        print(gg)
    }

    ## export
    return(invisible(list(data = CI,
                          plot = gg)))
}

## * .calcShortestCI
.calcShortestCI <- function(value, sigma, value.sign, threshold, alpha,
                            vec.init,
                            optimizer = "optim", method.optim = "L-BFGS-B", ...){

    ## ** find theta1 and theta2
    Q <- function(theta){ 2 - stats::pnorm(threshold + theta, mean = 0, sd = sigma) - stats::pnorm(threshold - theta, mean = 0, sd = sigma)	}
    d_Q <- function(theta){ - stats::dnorm(threshold + theta, mean = 0, sd = sigma) + stats::dnorm(threshold - theta, mean = 0, sd = sigma)	}

    fn1 <- function(theta){
        stats::pnorm(threshold + theta, mean = 0, sd = sigma) - stats::pnorm(threshold - theta, mean = 0, sd = sigma) - (1 - alpha) * Q(theta)
    }
    d_fn1 <- function(theta){
        stats::dnorm(threshold + theta, mean = 0, sd = sigma) + stats::dnorm(threshold - theta, mean = 0, sd = sigma) - (1 - alpha) * d_Q(theta)
    }
    fn2 <- function(theta){
        2*stats::pnorm(theta - threshold, mean = 0, sd = sigma) - 1 - (1 - alpha)* Q(theta)
    }
    d_fn2 <- function(theta){
        2*stats::dnorm(theta - threshold, mean = 0, sd = sigma) - (1 - alpha )* d_Q(theta)
    }

    if(optimizer == "optim"){
        iOut <- stats::optim(par = vec.init["theta1"],
                             fn = function(x){fn1(x)^2},
                             gr = function(x){2*fn1(x)*d_fn1(x)},
                             ## lower = 0,
                             method = method.optim)
        theta1.optim <- data.frame(solution = iOut$par,
                                   objective = iOut$value,
                                   iteration = iOut$counts["function"],
                                   convergence = iOut$convergence,
                                   message = iOut$message,
                                   stringsAsFactors = FALSE)
        theta1 <- as.double(iOut$par)

        ## curve(fn2, 0,10)
        ## method.optim <- "Nelder-Mead"
        iOut <- stats::optim(par = vec.init["theta2"],
                             fn = function(x){fn2(x)^2},
                             gr = function(x){2*fn2(x)*d_fn2(x)},
                             ## lower = theta1,
                             method = method.optim)
        theta2.optim <- data.frame(solution = iOut$par,
                                   objective = iOut$value,
                                   iteration = iOut$counts["function"],
                                   convergence = iOut$convergence,
                                   message = iOut$message,
                                   stringsAsFactors = FALSE)
        theta2 <- as.double(iOut$par)

    }else if(optimizer == "uniroot"){

        iOut <- stats::uniroot(f = fn1, lower = 0, upper = threshold + 5*qnorm(1 - alpha, mean = 0, sd = sigma))
        theta1.optim <- data.frame(solution = iOut$root,
                                   objective = iOut$f.root,
                                   iteration = iOut$iter,
                                   convergence = NA,
                                   message = "NA",
                                   stringsAsFactors = FALSE)
        theta1 <- iOut$root
        
        iOut <- stats::uniroot(f = fn2, lower = theta1, upper = threshold + 5*qnorm(1 - alpha/2, mean = 0, sd = sigma))
        theta2.optim <- data.frame(solution = iOut$root,
                                   objective = iOut$f.root,
                                   iteration = iOut$iter,
                                   convergence = NA,
                                   message = "NA",
                                   stringsAsFactors = FALSE)
        theta2 <- iOut$root
    }
    ## ** Identify acceptance region (optional, formula 2)
    if (value < theta1) { ## |theta| in [0,theta1]

        lower.left <- value - stats::qnorm( 1 - alpha/2 *  Q(value), mean = 0, sd = sigma)
        if(lower.left>=-threshold){
            lower.left <- NA
            upper.left <- NA
        }else{
            upper.left <- -threshold
        }
        lower.right <- threshold
        upper.right <- value + stats::qnorm( 1 - alpha/2 *  Q(value), mean = 0, sd = sigma)
        
    }else if (value < theta2) { ## |theta| in [theta1,theta2]
        
        lower.left <- NA
        upper.left <- NA
        lower.right <- threshold
        upper.right <- value + stats::qnorm( stats::pnorm( threshold - value, mean = 0, sd = sigma) + (1 - alpha) *  Q(value),
                                            mean = 0, sd = sigma)
        
    } else {  ## |theta| in [theta2,Inf]
        
        lower.left <- NA
        upper.left <- NA
        lower.right <- value - stats::qnorm( .5 * ( 1 + (1 - alpha) *  Q(value) ), mean = 0, sd = sigma)
        upper.right <- value + stats::qnorm( .5 * ( 1 + (1 - alpha) *  Q(value) ), mean = 0, sd = sigma)
        
    }

    ## reconstruct
    if (value.sign %in% c(0,1)){
        AR <- c(lower.left = lower.left, lower.right = lower.right, upper.left = upper.left, upper.right = upper.right)
    }else if (value.sign == -1){
        AR <- c(lower.left = -upper.right, lower.right = -upper.left, upper.left = -lower.right, upper.right = -lower.left)
    }

    ## ** Identify CI
    
    ## *** compute x1 and x2 according to formula (6-7)
    x1 <- threshold + 2*theta1
    x2 <- 2*theta2 - threshold

    ## *** lower CI, i.e. solve equation (5)
    if (threshold < value && value < x1) { ## x in [threshold, x1]
        f.lower <-  function(theta){
            2 * (1 - stats::pnorm( value - theta, mean = 0, sd = sigma)) - alpha * Q(theta)
        }
        d_f.lower <-  function(theta){
            2 * stats::dnorm( value - theta, mean = 0, sd = sigma) - alpha * d_Q(theta)
        }
    }else if (value < x2) { ## x in [x1, x2]
        f.lower <- function(theta){ 
            stats::pnorm( value - theta, mean = 0, sd = sigma) - stats::pnorm( threshold - theta, mean = 0, sd = sigma) - (1 - alpha) *  Q(theta)
        }
        d_f.lower <- function(theta){ 
            - stats::dnorm( value - theta, mean = 0, sd = sigma) + stats::dnorm( threshold - theta, mean = 0, sd = sigma) - (1 - alpha) *  d_Q(theta)
        }
    }else { ## x in [x2, Inf])
        f.lower <- function(theta){
            2 * stats::pnorm( value - theta, mean = 0, sd = sigma) - 1 - (1 - alpha) *  Q(theta)
        }
        d_f.lower <- function(theta){
            - 2 * stats::dnorm( value - theta, mean = 0, sd = sigma) - (1 - alpha) *  d_Q(theta)
        }
    }
        
    ## *** upper CI    
    f.upper <- function(theta){
        2 * stats::pnorm(theta - value, mean = 0, sd = sigma) - 1 - (1 - alpha) *  Q(theta)
    }
    d_f.upper <- function(theta){
        2 * stats::dnorm(theta - value, mean = 0, sd = sigma) - (1 - alpha) *  d_Q(theta)
    }

    if(optimizer == "optim"){
        iOut <- stats::optim(par = vec.init["lower"],
                             fn = function(x){f.lower(x)^2},
                             gr = function(x){2*f.lower(x)*d_f.lower(x)},
                             upper = value,
                             method = method.optim)
        lower.optim <- data.frame(solution = iOut$par,
                                  objective = iOut$value,
                                  iteration = iOut$counts["function"],
                                  convergence = iOut$convergence,
                                  message = iOut$message,
                                  stringsAsFactors = FALSE)
        lower <- as.double(iOut$par)

        iOut <- stats::optim(par = vec.init["upper"],
                             fn = function(x){f.upper(x)^2},
                             gr = function(x){2*f.upper(x)*d_f.upper(x)},
                             lower = max(lower,value),
                             method = method.optim)
        upper.optim <- data.frame(solution = iOut$par,
                                  objective = iOut$value,
                                  iteration = iOut$counts["function"],
                                  convergence = iOut$convergence,
                                  message = iOut$message,
                                  stringsAsFactors = FALSE)
        upper <- as.double(iOut$par)
        ## print(upper.optim$message)
        ## Q(vec.init["lower"])^2
    }else if(optimizer == "uniroot"){

        iOut <- stats::uniroot(f = f.lower, lower = -theta1, upper = value)
        lower.optim <- data.frame(solution = iOut$root,
                                  objective = iOut$f.root,
                                  iteration = iOut$iter,
                                  convergence = NA,
                                  message = "NA",
                                  stringsAsFactors = FALSE)
        lower <- as.double(iOut$root)
        
        iOut <- stats::uniroot(f = f.upper, lower = theta2, upper = threshold + 5*qnorm(1 - alpha/2, mean = 0, sd = sigma))
        upper.optim <- data.frame(solution = iOut$root,
                                  objective = iOut$f.root,
                                  iteration = iOut$iter,
                                  convergence = NA,
                                  message = "NA",
                                  stringsAsFactors = FALSE)
        upper <- as.double(iOut$root)
        ## uses Chebyshev's_inequality
        ## P[|X-mu|> 5*sd < 1/5^2]

    }

    ## M.print <- rbind(theta1 = unlist(theta1.optim[1:2]),
                     ## theta2 = unlist(theta2.optim[1:2]),
                     ## lower.optim = unlist(lower.optim[1:2]),
                     ## upper.optim = unlist(upper.optim[1:2])
                     ## )
    ## if(any(M.print[,"value"]>1e-7)){
        ## print(M.print)
    ## }
    
    ## *** reconstruct
    if (value.sign %in% c(0,1)){
        CI <- c(lower = lower, upper = upper)        
    }else if (value.sign == -1){
        CI <- c(lower = -upper, upper = -lower)
    }
    
    ## ** Export

    return(list(AR = AR,
                CI = CI,
                optim = rbind(cbind(param = "theta1", theta1.optim),
                              cbind(param = "theta2", theta2.optim),
                              cbind(param = "lower", lower.optim),
                              cbind(param = "upper", upper.optim))
                ))
}


######################################################################
### conditionalCI.R ends here

