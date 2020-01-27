### ggResPlot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 21 2020 (10:38) 
## Version: 
## Last-Updated: jan 21 2020 (11:13) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ggResPlot (documentation)
#' @title Residual plot using ggplot2
#' @description Display the residuals as a function of a variable
#' @name ggResPlot
#'
#' @param object a \code{lvm} object.
#' @param x [character] the name of an observed variable.
#' @param plot [logical] should the plot be displayed?.
#' @param geom [character] the type of plot to make:
#' can be \code{"point"} for a scatterplot or \code{"boxplot"} for boxplots.
#' @param ... arguments to be passed to lower level functions
#' @export
`ggResPlot` <-
  function(object, x, plot, ...) UseMethod("ggResPlot")

## * ggResPlot.lvmfit (code)
#' @rdname ggResPlot
#' @export
ggResPlot.lvmfit <- function(object, x, plot = TRUE, geom = "point", ...){

    ## check arguments
    if(x %in% manifest(object) == FALSE){
        stop("Argument \'x\' must correspond to a manifest variable in the lvmfit object \n")
    }
    if(any("XXindexXX" %in% manifest(object) == TRUE)){
        stop("\"XXindexXX\" is a named used internally and therefore should not be used to name a variable \n")
    }
    geom <- match.arg(geom, c("point","boxplot"))

    ## extract residuals
    df <- data.frame(object$data$model.frame[[x]],
                     stats::residuals(object))
    names(df)[1] <- x
    df$XXindexXX <- 1:NROW(df)

    ## reshape
    dtL <- data.table::melt(df, measure = setdiff(names(df),c("XXindexXX",x)) , id.var = c("XXindexXX",x),
                            variable.name = "endogenous", value.name = "residual")
    if(geom == "boxplot"){
        dtL[[x]] <- as.factor(dtL[[x]])
    }
    ## build plot
    gg <- ggplot2::ggplot(dtL, ggplot2::aes_string(x = x, y = "residual"))
    if(geom == "point"){
        gg <- gg + ggplot2::geom_point() + ggplot2::geom_smooth()
    }else if(geom == "boxplot"){
        gg <- gg + ggplot2::geom_boxplot()
    }
    gg <- gg + ggplot2::facet_wrap(~endogenous)

    ## display
    if(plot){
        print(gg)
    }

    ## export
    return(invisible(list(plot = gg,
                          data = dtL)))
}

######################################################################
### ggResPlot.R ends here
