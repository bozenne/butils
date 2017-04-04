# {{{ partialModel
# {{{ doc
#' @title Create a partial model relative to two variables
#' @description Create a partial model relative to two variables
#' @name partialModel
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

#' @rdname partialModel
#' @export
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
  model.matrix(terms(object), data = nlme::getData(object), ...)


model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = nlme::getData(object), ...)


terms.gls <- function(object, ...)
  terms(model.frame(object),...)
# }}}

# {{{ missing functions in lme 
model.matrix.lme <- function(object, ...)
  model.matrix(terms(object), data = nlme::getData(object), ...)


model.frame.lme <- function(object, ...)
  model.frame(formula(object), data = nlme::getData(object), ...)


terms.lme <- function(object, ...)
  terms(model.frame(object),...)
# }}}