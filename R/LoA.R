### LoA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2019 (09:46) 
## Version: 
## Last-Updated: apr  5 2019 (13:39) 
##           By: Brice Ozenne
##     Update #: 92
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * LOA - documentation
##' @title Limits of Agreement
##' @description Computes the limits of agreement between two measurements.
##' @name LoA
##'
##' @param x [numeric vector] vector containing the first measurement, each value corresponding to an individual.
##' @param y [numeric vector] vector containing the second measurement, each value corresponding to an individual (same as in \code{x}).
##' @param conf.level [numeric 0-1] confidence level of the confidence intervals.
##' @param LoA.multiplier [numeric 0-1] confidence level for the limits of agreements.
##'
##' @details The bias is estimated as the average difference between y and x, i.e. bias = E[y-x].
##' The upper/lower limit of the limit of agreement are estimated as the bias plus/minus the LoA.multiplier times the standard deviation of y-x,
##' i.e. E[y-x] +/- LoA.multiplier * sqrt(Var[y-x]).
##' 
##' The standard errors, confidence intervals, and p-values are computed for the limit of the limit of agreement are computed using a delta method.
##'
##' 
##' 
##' 
    
## * LoA - example
##' @rdname LoA
##' @examples
##'
##' set.seed(10)
##' X <- rnorm(100)
##' Y <- X + rnorm(100)
##'
##' e.loa <- LoA(x = X, y = Y)
##' e.loa
##' autoplot(e.loa)
##' 
##' if(require(blandr)){ ## same results as blandr.statistics
##'  e.blandr <- blandr.statistics(method1 = Y, method2 = X)
##'  e.loa$LoA$estimate - c(e.blandr$bias, e.blandr$lowerLOA, e.blandr$upperLOA)
##'  e.loa$LoA$lower - c(e.blandr$biasLowerCI, e.blandr$lowerLOA_lowerCI, e.blandr$upperLOA_lowerCI)
##'  e.loa$LoA$upper - c(e.blandr$biasUpperCI, e.blandr$lowerLOA_upperCI, e.blandr$upperLOA_upperCI)
##' }

## * LoA - code
##' @rdname LoA
##' @export
LoA <- function(x, y, conf.level = 0.95, LoA.multiplier = qnorm(0.975)){
        
    ## ** check arguments
    if(length(x) != length(y)){
        stop("Argument \'x\' and \'y\' must have same lenght \n")
    }

    ## ** normalize arguments
    alpha <- 1 - conf.level
    n <- length(x)
    df <- n-1
    critical.quantile <- qt(1-alpha/2, df = df)
    
    data <- data.frame(x = x,
                       y = y,
                       difference = y - x,
                       average = (y+x)/2)

    data.LoA <- data.frame(matrix(NA, nrow = 3, ncol = 5,
                                  dimnames = list(c("bias","lowerLoA","upperLoA"),
                                                  c("estimate","se","lower","upper","p.value"))))
    
    ## ** compute quantities of interest
    mean.diff <- mean(data$difference)
    
    var.diff <- stats::var(data$difference)
    ## var(sigma^2) = 2*sigma^4/(n-1)
    var.var.diff <- 2*var.diff^2/(n-1)
    
    sd.diff <- sqrt(var.diff)
    ## var(sigma) = var(sqrt(sigma^2)) = var(sigma^2)/(4*sigma^2) by a delta method
    var.sd.diff <- var.var.diff / (4*var.diff)
    
    se.diff <- sd.diff/sqrt(n)

    ## var(mu + q * sigma) = var(mu) + q^2 * var(sigma) + 2 * 0 (no covariance between mu and sigma)
    var.LoA <- (se.diff)^2 + LoA.multiplier^2 * var.sd.diff
    se.LoA <- sqrt(var.LoA)
    
    ## bias
    data.LoA["bias","estimate"] <- mean.diff
    data.LoA["bias","se"] <- se.diff

    ## lower LoA
    data.LoA["lowerLoA","estimate"] <- mean.diff - LoA.multiplier * sd.diff
    data.LoA["lowerLoA","se"] <- se.LoA

    ## upper LoA
    data.LoA["upperLoA","estimate"] <- mean.diff + LoA.multiplier * sd.diff
    data.LoA["upperLoA","se"] <- se.LoA

    ## p-value and CI
    data.LoA[,"lower"] <- data.LoA[,"estimate"] - critical.quantile * data.LoA[,"se"]
    data.LoA[,"upper"] <- data.LoA[,"estimate"] + critical.quantile * data.LoA[,"se"]
    data.LoA[,"p.value"] <- 2*(1-pt(abs(data.LoA[,"estimate"]/data.LoA[,"se"]), df = df))

    ## ** export
    out <- list(data = data,
                LoA = data.LoA)
    class(out) <- "LoA"
    return(out)
}

## * print.LoA
print.LoA <- function(x, digits = 2, eps = 1e-4, print = TRUE, ...){

    out <- x$LoA
    out$estimate <- round(out$estimate, digits = digits)
    out$se <- round(out$se, digits = digits)
    out$lower <- round(out$lower, digits = digits)
    out$upper <- round(out$upper, digits = digits)
    out$p.value <- format.pval(out$p.value, digits = digits, eps = eps)
    out$lower <- paste0("[",out$lower, ";",out$upper,"]")
    out$upper <- NULL
    names(out)[names(out)=="lower"] <- "CI"
    
    if(print == TRUE){
        print(out)
    }
    
    return(invisible(out))
}

## * autoplot.LoA
##' @title Bland and Altman Plot with Limits of Agreement
##' @description Display a Bland and Altman plot with limits of agreement.
##'
##' @param object An object of class \code{LoA}, output of the \code{LoA} function.
##' @param display.ci [logical] should the confidence intervals be displayed?
##' @param plot [logical] should the plot be printed?
##' @param line.size [numeric >0] the width of the horizontal lines.
##' @param alpha [numeric 0-1] transparency parameter for the confidence intervals.
##' @param name.legend [character] character string displayed in the caption.
##' @param colors [character of length 3] colors used to display the bias, lower, and upper LoA.
##' @param labels [character of length 3] text shown in the caption for the bias, lower, and upper LoA.
##' @param linetype [numeric of length 3] type of lines used to display the bias, lower, and upper LoA.
##' @param ...  not used, for compatibility with the generic method.
##' 
##' @export
autoplot.LoA <- function(object, display.ci = TRUE, plot = TRUE,
                         line.size = 1.5, alpha = 0.25, name.legend = "",
                         labels = c("bias","LoA (lower)","LoA (upper)"),
                         colors = c("red","blue","forestgreen"), linetype = c(1,2,2),
                         ...){


    if(display.ci){
        keep.type <- c("estimate","lower","upper")
    }else{
        keep.type <- c("estimate")
    }

    object$LoA$type <- rownames(object$LoA)

    ## display
    gg <- ggplot2::ggplot()
    if(display.ci){
        x.range <- range(object$data$average)
        gg <- gg + ggplot2::geom_rect(data = object$LoA, aes_string(ymin = "lower", ymax = "upper", fill = "type"), alpha = alpha, xmin = x.range[1], xmax = x.range[2])
    }

    gg <- gg + ggplot2::geom_point(data = object$data, mapping = ggplot2::aes_string(x = "average", y = "difference"))

    gg <- gg + ggplot2::geom_hline(data = data.frame(type = object$LoA$type, estimate = object$LoA$estimate),
                          ggplot2::aes_string(yintercept = "estimate", linetype = "type", color = "type"), size = line.size)

    gg <- gg + ggplot2::scale_colour_manual(name = name.legend,
                                   labels = labels,
                                   values = c("bias" = colors[1], 
                                              "lowerLoA" = colors[2],
                                              "upperLoA" = colors[3]))
    gg <- gg + ggplot2::scale_fill_manual(name = name.legend,
                                 labels = labels,
                                 values = c("bias" = colors[1], 
                                            "lowerLoA" = colors[2],
                                            "upperLoA" = colors[3]))
    
    gg <- gg + ggplot2::scale_linetype_manual(name = name.legend,
                                     labels = labels,
                                     values = c("bias" = linetype[1], 
                                                "lowerLoA" = linetype[2],
                                                "upperLoA" = linetype[3]))
    
    if(plot == TRUE){
        print(gg)
    }

    ## export
    return(invisible(gg))
}

######################################################################
### LoA.R ends here
