### FCT_plotInteraction.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr  4 2017 (15:45) 
## Version: 
## last-updated: maj 26 2017 (14:28) 
##           By: Brice Ozenne
##     Update #: 160
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc
#' @title Display an first order interaction for categorical variables
#' @description Display an first order interaction for categorical variables
#'
#' @param x a linear model
#' @param alpha.ggplot transparency parameter for the points or the bands
#' @param ... additional argument to be passed to \code{partialModel}
#'
#' @examples
#' library(lava)
#' set.seed(10)
#' m.lvm <- lvm(Y~2*Age+4*gender+gene+time)
#' categorical(m.lvm, labels = c("M","F")) <- ~gender
#' categorical(m.lvm, K = 10) <- ~gene
#' d <- sim(n = 1e2, m.lvm)
#'
#' ## linear model
#' m <- lm(Y~Age+gender, data = d)
#' pres1 <- calcPartialResiduals(m, var = "Age")
#' plot(pres1)
#' pres2 <- calcPartialResiduals(m, var = c("Age","gender"))
#' plot(pres2)
#' pres3 <- calcPartialResiduals(m, var = c("Age","gender"), interval = "prediction")
#' plot(pres3)
#' pres4 <- calcPartialResiduals(m, var = "gender")
#' plot(pres4)
#'
#' ## linear mixed model
#' if(require(lme4) && require(merTools) && require(AICcmodavg)){
#' 
#' mm <- lmer(Y~Age+gender+(1|gene), data = d)
#' pres1 <- calcPartialResiduals(mm, var = "Age")
#' plot(pres1)
#' pres2 <- calcPartialResiduals(mm, var = c("Age","gender"))
#' plot(pres2)
#'
#' # using external function
#' pres3 <- calcPartialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_merTools)
#' plot(pres3) 
#' pres4 <- calcPartialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_AICcmodavg)
#' plot(pres4)
#' 
#' }
# }}}


# {{{ plot.partialResiduals
#' @rdname plotInteraction
#' @export
plot.partialResiduals <- function(x, alpha.ggplot = 0.25, ...){    

    if(length(x$var)==1){
        var1 <- x$var[1]
        var2 <- NULL
    }else if(length(x$var)==2){
        var1 <- x$var[1]
        var2 <- x$var[2]
    }else if(length(x$var)>2){
        stop("the partial residuals must be computed only regarding 1 or 2 variables \n")
    }

    #### sort dataset
    x$data <- copy(x$data)
    x$partialFit <- cbind(x$partialFit, type = paste0(100*x$level,"% ",x$interval," interval"))
    if(is.null(var2)){
        setkeyv(x$data, var1)
        setkeyv(x$partialFit, var1)
    }else{
        setkeyv(x$data, c(var1,var2))
        setkeyv(x$partialFit, c(var1,var2))
    }
    
    #### display ####
    gg <- ggplot()
    if(is.null(var2) && !is.numeric(x$partialFit[[var1]])){
        gg <- gg + geom_point(data = x$data,
                              aes_string(x = var1, y = "pResiduals", color = var2), alpha = alpha.ggplot)
        gg <- gg + geom_point(data = x$partialFit,
                              aes_string(x = var1, y = "fit"), shape = 4, size = 4, color = "blue") 
        gg <- gg + geom_errorbar(data = x$partialFit,
                                 aes_string(x = var1, linetype = "type",  ymin = "fit.lower", ymax = "fit.upper"), color = "blue")
        gg <- gg + scale_linetype_manual("", values = 1)
        
    }else{
        gg <- gg + geom_point(data = x$data,
                              aes_string(x = var1, y = "pResiduals", color = var2))
        gg <- gg + geom_line(data = x$partialFit,
                             aes_string(x = var1, y = "fit", group = var2, color = var2))
        gg <- gg + geom_ribbon(data = x$partialFit,
                               aes_string(x = var1, fill = "type", ymin = "fit.lower", ymax = "fit.upper", group = var2, color = var2), alpha = alpha.ggplot)
        gg <- gg + scale_fill_manual("",values = "grey")
    }
    gg <- gg + ylab(paste0("Partial residuals \n",x$name.Y," | ",paste0(x$var,collapse = ", ")))
    print(gg)
     
    #### export ####
    return(invisible(gg))
    
}
# }}}



#----------------------------------------------------------------------
### FCT_plotInteraction.R ends here
