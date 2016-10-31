#' @title Box Cox transformation for linear models
#' @description Ouptut the linear model fitted on the Box-Cox transformed outcome
#' 
#' @param formula same as lm
#' @param data same as lm
#' @param trace should the execution of the function be traced
#' @param ... other arguments to be passed to lm
#'
#' @examples 
#' if(require(lava)){
#' df <- sim(lvm(Y ~ X1 + X2),1e2)
#' df$Y <- df$Y-min(df$Y)+1
#' lm(Y~X1+X2, data = df)
#' lmBC(Y~X1+X2, data = df)
#' #' 
#' }
#' 
#' @export
lmBC <- function(formula, data, trace = TRUE, ...){
  resBC <- MASS::boxcox(object = formula, data = data, ..., plot = FALSE)
  lambda <- resBC$x[which.max(resBC$y)]
  
  name.X <- all.vars(delete.response(terms(formula)))
  name.Y <- setdiff(all.vars(formula), name.X)
  if(trace){cat("boxcox transform with lambda = ",lambda,"\n")}
  if(lambda == 0){
    data[[name.Y]] <- log(data[[name.Y]])
  }else{
    data[[name.Y]] <- ((data[[name.Y]])^lambda-1) / lambda
  }
  lm.fit <- lm(formula = formula, data = data, ...)
  lm.fit$boxcox <- lambda
  return(lm.fit)
}



