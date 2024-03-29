### ggHeatmap.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  4 2018 (15:25) 
## Version: 
## Last-Updated: maj  4 2023 (13:09) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ggHeatmap (documentation)
#' @title Display pairwise correlation
#' 
#' @description Display a correlation matrix (or heatmap) using ggplot2
#' 
#' @param data the dataset containing the correlation values.
#' Either a matrix where each element is a correlation value or a data.frame/data.table containing in columns the correlation coefficients.
#' @param name.x name of the column containing the label of the first outcome involved in the correlation.
#' @param name.y name of the column containing the label of the second outcome involved in the correlation.
#' @param name.fill name of the column containing the value of the correlation.
#' @param plot should the graph be displayed.
#' @param add.text additional information to be displayed above the squares representing the correlation. Must name a column in data.
#' @param round if not \code{NULL}, the number of digit used to round add.text using signif. Can also be \code{p.value} to indicates with the usual convention the significance level.
#' @param title the title of the plot.
#' @param xlab the name of the x axis.
#' @param ylab the name of the y axis.
#' @param legend_title the name of the legend.
#' @param na.value the color used to display missing values (NA).
#' @param col_low the color used to diplay low correlation values,
#' @param col_midpoint the color used to diplay intermediate correlation values.
#' @param col_high the color used to diplay high correlation values.
#' @param midpoint the value corresponding to intermediate correlation value.
#' @param limits the minimum and maximum correlation values used to build the color panel.
#' @param text.size the size of the text in the plot.
#' @param angle.x the inclination of the x labels.
#' 
#' @details 
#' data must be coercible to data.table. \cr
#' 
#' @return a list containing:
#' \itemize{
#'   \item plot: a ggplot object.
#'   \item data: a data.table object used to create the plot.
#' }
#' @examples
#' 
#' fill <- runif(16)
#' df <- cbind(expand.grid(x = 1:4, y = 1:4),fill = fill)
#' ggHeatmap(df, name.x = "x", name.y = "y", name.fill = "fill")
#' 
#' M <- matrix(fill,4,4)
#' M[upper.tri(M)] <- NA
#' ggHeatmap(M)
#'
#' rownames(M) <- paste0("R",1:4)
#' ggHeatmap(M, limits = c(0,0.5), add.text = "fill", round = 2)
#' 
#' @keywords function correlation display

## * ggHeatmap (code)
#' @rdname ggHeatmap
#' @export
ggHeatmap <- function(data, name.x, name.y, name.fill, add.text, plot = TRUE, round = NULL,
                      title = "", xlab = "", ylab = "",  legend_title = "correlation",
                      na.value = "grey50", col_low = "blue", col_midpoint = "white", col_high = "red", midpoint = 0, limits = NULL, 
                      text.size = 15, angle.x = 90){
  XXIdXX <- NULL # For CRAN check
  
    if(is.matrix(data)){ # convert matrix to data.table
        
        if(!is.null(rownames(data))){
            row.names <- rownames(data)
        }else{
            row.names <- paste0("V",1:NROW(data))
        }
        level.x <- paste0("V",1:NROW(data))
        if(!is.null(colnames(data))){
            col.names <- colnames(data)
        }else{
            col.names <- paste0("V",1:NCOL(data))
        }
        level.y <- paste0("V",1:NCOL(data))
        if(missing(name.x)){
            name.x <- "x"
        }
        if(missing(name.y)){
            name.y <- "y"
        }
        if(missing(name.fill)){
            name.fill <- "fill"
        }

        data <- data.table::melt(data.table(XXIdXX = level.x, data), id.vars = "XXIdXX")
        data[,XXIdXX := factor(XXIdXX,  levels = level.x, labels = row.names)]
        #data[,variable := factor(variable,  levels = level.y, labels = col.names)]
        setnames(data, old = names(data),new = c(name.x,name.y,name.fill))

    }else{

        
        if(missing(name.x) && missing(name.y)){
            if(NROW(data)==NCOL(data)){
                stop("Arguments \'name.x\' and \'name.y\' are missing \n",
                     "Consider applying as.matrix() to the argument \'data\' \n")
            }else{
                stop("Arguments \'name.x\' and \'name.y\' are missing but are required when argument \'data\' is a \"",class(data)[1],"\"\n")
            }
        }
        
        if(!is.data.table(data)){
            data <- as.data.table(data)
        }else{
            data <- copy(data)
        }
    }

    ## ** prepare
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
  }else{
      data[[name.fill]][data[[name.fill]]<limits[1]] <- limits[1]
      data[[name.fill]][data[[name.fill]]>limits[2]] <- limits[2]
  }

  ## ** plot
  gg <- ggplot2::ggplot(data, aes_string(x = name.x, y = name.y)) 
  gg <- gg + ggplot2::geom_tile(aes_string(fill = name.fill))
  gg <- gg + ggplot2::ggtitle(title)
  gg <- gg + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) 
  gg <- gg + ggplot2::scale_y_discrete(limits = rev(levels(data[[name.y]])))
  gg <- gg + ggplot2::scale_fill_gradient2(low = col_low,
                                           mid = col_midpoint,
                                           high = col_high, 
                                           name = legend_title, 
                                           midpoint = midpoint,
                                           na.value = na.value,
                                           limits = limits)
  gg <- gg + ggplot2::theme(text = element_text(size=text.size),
                            axis.text.x = element_text(angle = angle.x, hjust = 1))
  gg <- gg + ggplot2::theme(legend.key.height = unit(0.1,"npc"),
                            legend.key.width = unit(0.08,"npc"))
  if(!missing(add.text)){
      gg <- gg + ggplot2::geom_text(aes_string(x = name.x, y = name.y, label = "add.text"))  
  }
    
  ## ** export
  if(plot == TRUE){
      print(gg)
  }
  
  return(invisible(list(plot = gg,
                        data = data)))
}






######################################################################
### ggHeatmap.R ends here
