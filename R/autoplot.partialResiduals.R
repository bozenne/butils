### FCT_plotInteraction.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr  4 2017 (15:45) 
## Version: 
## last-updated: feb 14 2020 (11:03) 
##           By: Brice Ozenne
##     Update #: 182
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - autoplot.partialResiduals
#' @title Display an first order interaction for categorical variables
#' @description Display an first order interaction for categorical variables
#'
#' @rdname autoplot.partialResiduals
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
#' d <- lava::sim(n = 1e2, m.lvm)
#'
#' ## linear model
#' m <- lm(Y~Age+gender, data = d)
#' pres1 <- partialResiduals(m, var = "Age")
#' autoplot(pres1)
#' pres2 <- partialResiduals(m, var = c("Age","gender"))
#' autoplot(pres2)
#' pres3 <- partialResiduals(m, var = c("Age","gender"), interval = "prediction")
#' autoplot(pres3)
#' pres4 <- partialResiduals(m, var = "gender")
#' autoplot(pres4)
#'
#' ## linear mixed model
#' if(require(lme4) && require(merTools) && require(AICcmodavg)){
#' 
#' mm <- lmer(Y~Age+gender+(1|gene), data = d)
#' pres1 <- partialResiduals(mm, var = "Age")
#' autoplot(pres1)
#' pres2 <- partialResiduals(mm, var = c("Age","gender"))
#' autoplot(pres2)
#'
#' # using external function
#' pres3 <- partialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_merTools)
#' autoplot(pres3) 
#' pres4 <- partialResiduals(mm, var = c("Age","gender"), FUN.predict = predict_AICcmodavg)
#' autoplot(pres4)
#' }


## * autoplot.partialResiduals
#' @rdname autoplot.partialResiduals
#' @method autoplot partialResiduals
#' @export
autoplot.partialResiduals <- function(object, alpha.ggplot = 0.25, ...){    
        
    name.Y <- attr(object,"name.Y")
    object.var <- attr(object,"var")
    if(length(object.var)==1){
        var1 <- object.var[1]
        var2 <- NULL
    }else if(length(object.var)==2){
        var1 <- object.var[1]
        var2 <- object.var[2]
    }else{
        stop("the partial residuals must be computed only regarding 1 or 2 variables \n")
    }

    ## ** sort dataset
    data <- copy(object)
    partialFit <- cbind(attr(object,"partialFit"), type = paste0(100*attr(object,"level"),"% ",attr(object,"interval")," interval"))
    if(is.null(var2)){
        setkeyv(data, var1)
        setkeyv(partialFit, var1)
    }else{
        setkeyv(data, c(var1,var2))
        setkeyv(partialFit, c(var1,var2))
    }
    
### ** display
    gg <- ggplot2::ggplot()
    if(is.null(var2) && !is.numeric(partialFit[[var1]])){
        gg <- gg + ggplot2::geom_point(data = data,
                                       aes_string(x = var1, y = "pResiduals", color = var2),
                                       alpha = alpha.ggplot)
        gg <- gg + ggplot2::geom_point(data = partialFit,
                                       aes_string(x = var1, y = "fit"),
                                       shape = 4,
                                       size = 4,
                                       color = "blue") 
        gg <- gg + ggplot2::geom_errorbar(data = partialFit,
                                          aes_string(x = var1, linetype = "type",  ymin = "fit.lower", ymax = "fit.upper"),
                                          color = "blue")
        gg <- gg + ggplot2::scale_linetype_manual("", values = 1)
        
    }else{
        gg <- gg + ggplot2::geom_point(data = data,
                                       aes_string(x = var1, y = "pResiduals", color = var2))
        gg <- gg + ggplot2::geom_line(data = partialFit,
                                      aes_string(x = var1, y = "fit", group = var2, color = var2))
        gg <- gg + ggplot2::geom_ribbon(data = partialFit,
                                        aes_string(x = var1, fill = "type", ymin = "fit.lower", ymax = "fit.upper", group = var2, color = var2),
                                        alpha = alpha.ggplot)
        gg <- gg + ggplot2::scale_fill_manual("",values = "grey")
    }
    gg <- gg + ggplot2::ylab(paste0("Partial residuals \n",name.Y," | ",paste0(object.var,collapse = ", ")))
    gg <- gg + ggplot2::theme(legend.position = "bottom")
    print(gg)
     
## ** export
    return(invisible(gg))
    
}



#----------------------------------------------------------------------
### FCT_plotInteraction.R ends here
