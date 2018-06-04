## * ggCor (documentation)
#' @title Display correlation/covariance  matrix
#' @description Display correlation/covariance matrix
#' @name ggCor
#' 
#' @param data a matrix containing in each columns the observations for each variable.
#' @param type should the raw values, correlation or covariance be displayed?
#' @param legend_title argument to be passed to \code{\link{ggHeatmap}}.
#' @param plot should the graphic be plotted?
#' @param ... additionnal parameters to be passed to \code{\link{ggHeatmap}}.
#' 
#' @examples 
#' 
#' n <- 100
#' X <- cbind(X1=rnorm(n),X2=rnorm(n),X3=rnorm(n))
#' 
#' ggCor(X)

## * ggCor (code)
#' @rdname ggCor
#' @export
ggCor <- function(data, type = "correlation", legend_title = NULL, plot = TRUE, ...){
  
  data <- as.data.frame(data)
  
  match.arg(arg = type, choices = c("none","correlation","covariance"))
  if(type %in% c("correlation","covariance")){
    data2 <- cov(data)
    if(type == "correlation"){
      data2 <- cov2cor(data2)
    }
  }else{
    data2 <- data
  }
  
  if(is.null(legend_title)){
    if(type == "none"){
      legend_title <- "value"
    }else{
      legend_title <- paste0(type,"\n")
    }
  }
  
  data2 <- data.table::melt(data2)
  out <- ggHeatmap(data2, 
                   name.x = "Var1", name.y = "Var2", name.fill = "value", 
                   legend_title = legend_title, ...)
  if(plot){
    print(out)
  }
  
  return(invisible(list(plot = out,
                        data = data2)))
                    
}
                    
