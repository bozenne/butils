### FCT_plotInteraction.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr  4 2017 (15:45) 
## Version: 
## last-updated: sep 27 2021 (17:08) 
##           By: Brice Ozenne
##     Update #: 252
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
#' @name autoplot.partialResiduals
#' 
#' @param x a linear model
#' @param size.point [numeric, >0] size of the dots representing the observed data.
#' @param col.point [character vector] color of the dots representing the observed data.
#' @param shape.fit [integer, >0] Symbol used to represent the fitted value.
#' @param size.fit [numeric, >0] thickness of the regression line.
#' @param col.fit [character vector] color of the regression line.
#' @param size.ci [numeric, >0] thickness of the line representing the confidence interval.
#' @param alpha.ggplot [numeric, 0-1] transparency parameter for the confidence interval.
#' @param col.ci [character vector] color of the line/band representing the confidence interval.
#' @param plot [logical]should the plot be displayed?
#' @param ... ignored.
#'
#' @examples
#' library(lava)
#' set.seed(10)
#' m.lvm <- lvm(Y~2*Age+4*gender+gene+time)
#' categorical(m.lvm, labels = c("M","F")) <- ~gender
#' categorical(m.lvm, K = 10) <- ~gene
#' d <- lava::sim(n = 1e2, m.lvm)
#' d$gene <- as.character(d$gene)
#'
#' ## linear model
#' m <- lm(Y~Age+gender, data = d)
#' pres1 <- partialResiduals(m, var = "Age")
#' autoplot(pres1)
#' autoplot(pres1, col.point = "pink", col.ci = "orange", col.fit = "purple")
#' pres2 <- partialResiduals(m, var = c("Age","gender"))
#' autoplot(pres2)
#' pres3 <- partialResiduals(m, var = c("Age","gender"), interval = "prediction")
#' autoplot(pres3)
#' pres4 <- partialResiduals(m, var = "gender")
#' autoplot(pres4)
#' autoplot(pres4, col.point = "red", col.ci = "orange", col.fit = "purple")
#' m2 <- lm(Y~Age+gender+gene, data = d)
#' pres5 <- partialResiduals(m2, var = c("gender","gene"))
#' autoplot(pres5)
#' autoplot(pres5, dodge = 1, col.point = rep("black",10))
#'
#' ## linear mixed model
#' if(require(nlme)){
#' mm <- lme(Y~Age+gender, random= ~ 1|gene, data = d)
#' pres1 <- partialResiduals(mm, var = "Age")
#' autoplot(pres1)
#' }
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
#'
#' ## gam
#' if(require(mgcv)){
#' set.seed(2) ## simulate some data
#' dat <- gamSim(1,n=400,dist="normal",scale=2)
#' b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
#' b <- gam(y~s(x0)+x1+x2+x3,data=dat)
#' pres5 <- partialResiduals(b, var = "x0")
#' autoplot(pres5)
#' }


## * autoplot.partialResiduals
#' @rdname autoplot.partialResiduals
#' @method autoplot partialResiduals
#' @export
autoplot.partialResiduals <- function(object,
                                      size.point = 2, col.point = NULL, dodge = 0.1,
                                      size.fit = NULL, shape.fit = NULL, col.fit = NULL,
                                      size.ci = 0.25, alpha.ggplot = 0.25, col.ci = NULL,
                                      plot = TRUE, ...){    

    ## ** normalize user input
    if(length(list(...))>0){
        warning("Argument \"",paste(names(list(...)), collapse = "\" \""),"\" will be ignored. \n")
    }
    
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

    ## ** display
    gg <- ggplot2::ggplot()
    if(is.null(var2)){
        if(is.null(col.point)){col.point <- "black"}
        if(is.null(col.fit)){col.fit <- "blue"}
        if(is.null(col.ci)){col.ci <- "blue"}
        
        if(is.numeric(partialFit[[var1]])){ ## continuous outcome, i.e. scatterplot
            if(is.null(size.fit)){size.fit <- 1.25}
            if(is.null(shape.fit)){shape.fit <- 1}

            gg <- gg + ggplot2::geom_point(data = data,
                                           aes_string(x = var1, y = "pResiduals"), size = size.point, color = col.point)
            gg <- gg + ggplot2::geom_line(data = partialFit,
                                          aes_string(x = var1, y = "fit"), size = size.fit, linetype = shape.fit, color = col.fit)
            gg <- gg + ggplot2::geom_ribbon(data = partialFit,
                                            aes_string(x = var1, ymin = "fit.lower", ymax = "fit.upper"),
                                            alpha = alpha.ggplot, size = size.ci, fill = col.ci)
        
        }else{ ## binary outcome, i.e. boxplot
            if(is.null(size.fit)){size.fit <- 4}
            if(is.null(shape.fit)){shape.fit <- 4}
            gg <- gg + ggplot2::geom_point(data = data,
                                           aes_string(x = var1, y = "pResiduals"),
                                           alpha = alpha.ggplot, size = size.point, color = col.point)
            gg <- gg + ggplot2::geom_point(data = partialFit,
                                           aes_string(x = var1, y = "fit"),
                                           shape = shape.fit,
                                           color = col.fit, size = size.fit) 
            gg <- gg + ggplot2::geom_errorbar(data = partialFit,
                                              aes_string(x = var1, linetype = "type",  ymin = "fit.lower", ymax = "fit.upper"),
                                              color = col.ci, size = size.ci)
            gg <- gg + ggplot2::scale_linetype_manual("", values = 1)            
        }
    }else{
        
        if(is.numeric(partialFit[[var1]])){ ## continuous outcome, i.e. scatterplot

            if(is.null(size.fit)){size.fit <- 1.25}
            if(is.null(shape.fit)){shape.fit <- 1}
            if(is.null(col.ci)){col.ci <- "black"}
            if(!is.null(col.fit) && !identical(col.point,col.fit)){
                stop("For partial residuals relative to two covariates, argument \'col.point\' and \'col.fit\' should be identical or NULL. \n")
            }
            gg <- gg + ggplot2::geom_point(data = data,
                                           aes_string(x = var1, y = "pResiduals", color = var2), size = size.point)
            gg <- gg + ggplot2::geom_line(data = partialFit,
                                          aes_string(x = var1, y = "fit", group = var2, color = var2), size = size.fit, linetype = shape.fit)
            gg <- gg + ggplot2::geom_ribbon(data = partialFit,
                                            aes_string(x = var1, fill = "type", ymin = "fit.lower", ymax = "fit.upper", group = var2, color = var2),
                                            alpha = alpha.ggplot, size = size.ci)
            if(!is.null(col.point)){
                gg <- gg + ggplot2::scale_color_manual("",values = col.point)
            }
            gg <- gg + ggplot2::scale_fill_manual("",values = col.ci)
        }else{
            if(is.null(size.fit)){size.fit <- 4}
            if(is.null(shape.fit)){shape.fit <- 4}
            if(!is.null(col.fit) && !identical(col.point,col.fit)){
                stop("For partial residuals relative to two covariates, argument \'col.point\' and \'col.fit\' should be identical or NULL. \n")
            }
            if(!is.null(col.ci) && !identical(col.point,col.ci)){
                stop("For partial residuals relative to two covariates, argument \'col.point\' and \'col.ci\' should be identical or NULL. \n")
            }
            ## gg <- ggplot2::ggplot()
            gg <- gg + ggplot2::geom_point(data = data,
                                           aes_string(x = var1, y = "pResiduals", color = var2),
                                           alpha = alpha.ggplot, size = size.point, position = position_dodge(width = dodge))
            gg <- gg + ggplot2::geom_point(data = partialFit,
                                           aes_string(x = var1, y = "fit", color = var2),
                                           shape = shape.fit,
                                           size = size.fit, position = position_dodge(width = dodge)) 
            gg <- gg + ggplot2::geom_errorbar(data = partialFit,
                                              aes_string(x = var1, linetype = "type",  ymin = "fit.lower", ymax = "fit.upper", color = var2),
                                              size = size.ci, position = position_dodge(width = dodge))
            gg <- gg + ggplot2::scale_linetype_manual("", values = 1)            
            if(!is.null(col.point)){
                gg <- gg + ggplot2::scale_color_manual("",values = col.point)
            }
        }
    }
        
    gg <- gg + ggplot2::ylab(paste0("Partial residuals \n",name.Y," | ",paste0(object.var,collapse = ", ")))
    gg <- gg + ggplot2::theme(legend.position = "bottom")

    if(plot){
        print(gg)
    }
    
    ## ** export
    return(invisible(gg))
    
}



#----------------------------------------------------------------------
### FCT_plotInteraction.R ends here
