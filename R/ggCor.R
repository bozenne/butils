#' @title Display correlation/covariance  matrix
#' @description Display correlation/covariance matrix
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
#' 
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
  
  data2 <- melt(data2)
  out <- ggHeatmap(data2, 
                   name.x = "Var1", name.y = "Var2", name.fill = "value", 
                   legend_title = legend_title, ...)
  if(plot){
    print(out)
  }
  
  return(invisible(list(plot = out,
                        data = data2)))
                    
}
                    

#' @title Display pairwise correlation
#' 
#' @description Display a correlation matrix (or heatmap) using ggplot2
#' 
#' @param data the dataset containing the correlation values.
#' @param name.x name of the column containing the label of the first outcome involved in the correlation
#' @param name.y name of the column containing the label of the second outcome involved in the correlation
#' @param name.fill name of the column containing the value of the correlation
#' @param add.text additional information to be displayed above the squares representing the correlation. Must name a column in data.
#' @param round if not \code{NULL}, the number of digit used to round add.text using signif. Can also be \code{p.value} to indicates with the usual convention the significance level.
#' @param title the title of the plot
#' @param xlab the name of the x axis
#' @param ylab the name of the y axis
#' @param legend_title the name of the legend
#' @param na.value the color used to display missing values (NA).
#' @param col_low the color used to diplay low correlation values
#' @param col_midpoint the color used to diplay intermediate correlation values
#' @param col_high the color used to diplay high correlation values
#' @param midpoint the value corresponding to intermediate correlation value.
#' @param limits the minimum and maximum correlation values used to build the color panel
#' @param textSize the size of the text in the plot
#' @param angle.x the inclination of the x labels
#' 
#' @details 
#' data must be coercible to data.table. \cr
#' 
#' @return a ggplot object
#' 
#' @keywords function correlation display
#' @export
ggHeatmap <- function(data, name.x, name.y, name.fill, add.text, round = NULL,
                      title = "", xlab = "", ylab = "",  legend_title = "correlation",
                      na.value = "grey50", col_low = "blue", col_midpoint = "white", col_high = "red", midpoint = 0, limits = NULL, 
                      textSize = 15, angle.x = 90){
  
  if(!is.data.table(data)){
    data <- as.data.table(data)
  }else{
    data <- copy(data)
  }
  
  #### prepare
  if(!missing(add.text)){
    
    if(add.text %in% names(data) == FALSE){
      stop("ggHeatmap: variable ",add.text," not found \n",
           "add.text should be one of \"",paste(names(data), collapse = "\" \""),"\"\n")
    }
    
    if(!is.null(round)){
      if(is.numeric(round)){
        data[, add.text := round(.SD[[1]], digit = round), .SDcols = add.text]
      }else if(round == "p.value"){
        type <- findInterval(data[["p.value"]], vec = c(0.001,0.01,0.05,0.1))
        type <- factor(type, levels = 0:4, labels = c("***","**","*",".",""))
        data[, add.text := type]
      }else{
        stop("ggHeatmap: non valid value for argument \'round\' \n")
      }
    }else{
      data[, add.text := .SD[[1]], .SDcols = add.text]
    }
    
  }
  
  data[[name.x]] <- factor(data[[name.x]], levels = unique(data[[name.x]]))
  data[[name.y]] <- factor(data[[name.y]], levels = unique(data[[name.y]]))
  
  if(is.null(limits)){
    limits <- range(data[[name.fill]], na.rm = TRUE)
  }
  
  #### plot
  gg <- ggplot(data, aes_string(x = name.x, y = name.y, fill = name.fill)) + geom_tile()
  gg <- gg + ggtitle(title)
  gg <- gg + xlab(xlab) + ylab(ylab) +  scale_y_discrete(limits = rev(levels(data[[name.x]])))
  gg <- gg + scale_fill_gradient2(low = col_low, mid = col_midpoint, high = col_high, 
                                  name = legend_title, 
                                  midpoint = midpoint, na.value = na.value, limits = limits)
  gg <- gg + theme(text = element_text(size=textSize), axis.text.x = element_text(angle = angle.x, hjust = 1))
  gg <- gg + theme(legend.key.height=unit(0.1,"npc"),
                   legend.key.width=unit(0.08,"npc"))
  if(!missing(add.text)){
    gg <- gg +  geom_text(aes_string(fill = name.fill, label = "add.text"))  
  }
  
  return(gg)
}




