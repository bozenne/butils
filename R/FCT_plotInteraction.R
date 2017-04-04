### FCT_plotInteraction.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr  4 2017 (15:45) 
## Version: 
## last-updated: apr  4 2017 (18:10) 
##           By: Brice Ozenne
##     Update #: 104
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ plotInteraction
# {{{ doc
#' @title Display an first order interaction for categorical variables
#' @description Display an first order interaction for categorical variables
#'
#' @param object a linear model
#' @param var1 a continuous or categorical covariate
#' @param var2 a categorical covariate
#' @param method.predict a method to compute the predicted value from the model (fit) with their standard error (se.fit).
#' @param field.data the name of the data field in the call
#' @param alpha the alpha risk for the confidence interval
#' @param n.points the number of points use for the predictions along the x axis. Only relevant if var1 is numerical.
#' @param plot should the plot be displayed.
#' @param alpha.ggplot the transparency level for the display of the confidence intervals.
#' @param ... additional argument to be passed to \code{partialModel}
#'
#' @examples
#' library(lava)
#'
#' ## simulate interaction
#' m <- lvm(Y ~ X1 + gender + group + Interaction)
#' distribution(m, ~gender) <- binomial.lvm()
#' distribution(m, ~group) <- binomial.lvm()
#' constrain(m, Interaction ~ gender + group) <- function(x){x[,1]*x[,2]}
#' d <- sim(m, 1e2)
#' d$gender <- factor(d$gender, labels = letters[1:2])
#' d$group <- factor(d$group)
#' 
#' ## linear model
#' lmfit <- lm(Y ~ X1 + gender*group, data = d)
#'
#' plotInteraction(lmfit, var1 = "gender", var2 = "group")
#' plotInteraction(lmfit, var1 = "X1", var2 = "group")
#'
#' ## gls model
#' glsfit <- gls(Y ~ X1 + gender*group, data = d,
#'              weight = varIdent(form = ~1|group))
#'
#' plotInteraction(glsfit, var1 = "gender", var2 = "group")
#' plotInteraction(glsfit, var1 = "X1", var2 = "group")
#'
#' ## mixed model
#' lmefit <- lme(Y ~ X1 + gender*group, data = d,
#'               random = ~1|group)
#'
#' plotInteraction(lmefit, var1 = "gender", var2 = "group")
#' res <- plotInteraction(lmefit, var1 = "X1", var2 = "group")
#' res$model
#'
# }}}

plotInteraction <- function(object, var1, var2,
                            method.predict = NULL, field.data = NULL,
                            alpha = 0.05, n.points = 1000, plot = TRUE, alpha.ggplot = 0.25, ...){

    if(is.null(method.predict)){
        if(any(class(object) %in% c("gls","lme"))){
            if(require(AICcmodavg)){
                method.predict <- "predictSE"
            }else{
                stop("the package AICcmodavg needs to be installed \n")
            }
        }else{
            method.predict <- "predict"
        }        
    }
    if(is.null(field.data)){
        field.data <- "data"
    }

    
    #### 1- refit model while removing the effect of the other covariates ####
    object <- partialModel(object, var1 = var1, var2 = var2, ...)

    #### 2- form dataset for prediction ####
    
    ## get outcome name
    name.Y <- all.vars(formula(object)[[2]])

    ## get data
    data <- eval(object$call[[field.data]])

    ## form newdata
    if(is.numeric(data[[var1]])){
        seqVar1 <- seq(min(data[[var1]], na.rm = TRUE), max(data[[var1]], na.rm = TRUE), length.out = n.points)
    }else{
        seqVar1 <- unique(data[[var1]])
    }    
    if(is.numeric(data[[var2]])){
        stop("The variable corresponding to argument \'var2\' should be categorical \n")
    }
    
    seqVar2 <- unique(data[[var2]])
    newdata <- expand.grid(seqVar1, seqVar2)
    names(newdata) <- c(var1,var2)

    #### 3- perform predictions ####
    res.predict <- do.call(method.predict, args = list(object, newdata = newdata, se.fit = TRUE))
    newdata$fit <- res.predict$fit

    if("df" %in% names(res.predict)){
        qt <- c(lower = qt(alpha/2, df = res.predict$df),
                upper = qt(1-alpha/2, df = res.predict$df))
    }else{
        qt <- c(lower = qnorm(alpha/2),
                upper = qnorm(1-alpha/2))
    }
    newdata$fit <- res.predict$fit
    newdata$fit.lower <- res.predict$fit + qt["lower"]*res.predict$se.fit
    newdata$fit.upper <- res.predict$fit + qt["upper"]*res.predict$se.fit
    
    #### 4- display ####
    gg <- ggplot()
    gg <- gg + geom_point(data = data, aes_string(x = var1, y = name.Y, color = var2))
    gg <- gg + geom_line(data = newdata, aes_string(x = var1, y = "fit", group = var2, color = var2))
    gg <- gg + geom_ribbon(data = newdata, aes_string(x = var1, ymin = "fit.lower", ymax = "fit.upper", group = var2, color = var2), alpha = alpha.ggplot)
    if(plot){
        print(gg)
    }

    #### export ####
    ls.export <- list(plot = gg,
                      data =  newdata,
                      model = object)
    return(invisible(ls.export))
    
}
# }}}

# {{{ partialModel
# {{{ doc
#' @title Create a partial model relative to two variables
#' 
#' @param object the model to be reduced
#' @param var1 a covariate
#' @param var2 another covariate
#' @param method.coef the method to extract the coefficient from the model
#' @param field.formula the name of the formula field in the call
#' @param field.data the name of the data field in the call
#' 
#' @return the reduced model
#' 
#' @examples
#' library(lava)
#' library(nlme)
#' 
#' ## lm
#' m <- lvm(Y~X1+X2+X3+X4+X5)
#' categorical(m, K=3) <- ~X1
#' categorical(m, K=2, labels = letters[1:2]) <- ~X2
#' d <- sim(m, 1e2)
#'
#' lm.full <- lm(Y~X1*X2+X3+X4+X5, data = d)
#' coef(lm.full)
#' partialModel(lm.full, var1 = "X4", var2 = "X5")
#' partialModel(lm.full, var1 = "X1", var2 = "X5")
#' partialModel(lm.full, var1 = "X2", var2 = "X3")
#' partialModel(lm.full, var1 = "X1", var2 = "X2")
#'
#' gls.full <- gls(Y~X1*X2+X3+X4+X5, data = d,
#'                 weights = varIdent(form = ~ 1|X2))
#' coef(gls.full)
#' partialModel(gls.full, var1 = "X4", var2 = "X5")
#' partialModel(gls.full, var1 = "X1", var2 = "X5")
#' partialModel(gls.full, var1 = "X2", var2 = "X3")
#' partialModel(gls.full, var1 = "X1", var2 = "X2")
#'
#' lme.full <- lme(Y~X1*X3+X4+X5, data = d, random = ~ 1|X2)
#' fixef(lme.full)
#' partialModel(lme.full, var1 = "X4", var2 = "X5")
#' partialModel(lme.full, var1 = "X1", var2 = "X5")
#' partialModel(lme.full, var1 = "X2", var2 = "X3")
#' partialModel(lme.full, var1 = "X1", var2 = "X2")
#' 
# }}}

partialModel <- function(object, var1, var2,
                         method.coef = NULL, field.formula = NULL, field.data = NULL){

    ## normalize model sturcture
    if(is.null(method.coef)){
        if(any(class(object) %in% "lme")){
            method.coef <- "fixef"
        }else{
            method.coef <- "coef"
        }
        
    }
    if(is.null(field.formula)){
        if(any(class(object) %in% "gls")){
            field.formula <- "model"
        }else if(any(class(object) %in% "lme")){
            field.formula <- "fixed"
        }else{
            field.formula <- "formula"
        }
    }
    if(is.null(field.data)){
        field.data <- "data"
    }

    ## get formula
    formula.object <- formula(object)
    
    ## get outcome name
    name.Y <- all.vars(formula.object[[2]])
        
    ## get coef
    Beta <- do.call(method.coef, args = list(object))

    ## get data
    data <- eval(object$call[[field.data]])
    X <- model.matrix(formula.object, data)

    ## index remove columns corresponding to var1 and var2
    name.rm <- list()
    if(is.numeric(data[,var1])){
        name.rm[[1]] <- var1
    }else{
        name.rm[[1]] <- paste0(var1,unique(levels(data[,var1])))
    }
    if(is.numeric(data[[var2]])){
        name.rm[[2]] <- var2
    }else{
        name.rm[[2]] <- paste0(var2,unique(levels(data[,var2])))
    }

    interaction.rm <-  apply(expand.grid(name.rm[[1]], name.rm[[2]]), 1,
                             paste, collapse = ":")

    index.rm <- which(names(Beta) %in% c(name.rm[[1]], name.rm[[2]], interaction.rm))

    if(length(index.rm)>0){ ## create the partial model
        data[[name.Y]] <- data[[name.Y]] - X[,-index.rm,drop=FALSE] %*% Beta[-index.rm]

        if(any(interaction.rm %in% names(Beta))){
            txt.formula <- paste0(".~",var1,"*",var2)
        }else{
            txt.formula <- paste0(".~",var1,"+",var2)
        }

        newFormula <- update(formula.object, as.formula(txt.formula))
        ls.update <- list(object, newFormula, data = data)
        names(ls.update)[2] <- field.formula
        object <- do.call(update, ls.update)
        object <- do.call(update, list(object, data = data))        
    }

    ## export
    newBeta <- do.call(method.coef, args = list(object))
    newBeta <- newBeta[names(newBeta) %in% c(name.rm[[1]], name.rm[[2]], interaction.rm)]
    if( any(abs(newBeta-Beta[names(newBeta)])>1e-6) ){
        warning("The coefficients of the partial model seems to differ compare to the original model \n",
                "range of the difference: ",paste(range(newBeta-Beta[names(newBeta)]), collapse = " ; "),"\n")
    }
    return(object)

}
# }}}

# {{{ missing functions in gls 
model.matrix.gls <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)


model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)


terms.gls <- function(object, ...)
  terms(model.frame(object),...)
# }}}

# {{{ missing functions in lme 
model.matrix.lme <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)


model.frame.lme <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)


terms.lme <- function(object, ...)
  terms(model.frame(object),...)
# }}}



#----------------------------------------------------------------------
### FCT_plotInteraction.R ends here
