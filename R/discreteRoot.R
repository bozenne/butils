### discreteRoot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (13:39) 
## Version: 
## Last-Updated: nov 22 2017 (16:33) 
##           By: Brice Ozenne
##     Update #: 67
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Dichotimic search for monotone function
#' @description Find the root of a monotone function on a discrete grid of value using dichotomic search
#'
#' @param fn objective function to minimize in absolute value
#' @param grid a vector of values
#' @param increasing is the function fn increasing?
#' @param check should the program check that fn takes a different sign for the first vs. the last value of the grid?
#' @param tol The absolute convergence tolerance
#'
#' @examples
#'
#' ### find the position of a value in a vector
#' f <- function(x){abs(vec[x]-1)}
#' discreteRoot(function(x){x},grid = seq(-20,10,1))
#' 
#' ### find level of the confidence interval
#' library(nlme)
#' fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
#'                correlation = corAR1(form = ~ 1 | Mare))
#'
#' fctIC1 <- function(x){    
#'    IC.tempo <- intervals(fm1, level = 1-x)
#'    return( IC.tempo[["coef"]][1,"upper"])
#' }
#' fctIC2 <- function(x){    
#'    IC.tempo <- intervals(fm1, level = 1-x)
#'    return( IC.tempo[["coef"]][2,"upper"])
#' }
#' fctIC3 <- function(x){    
#'    IC.tempo <- intervals(fm1, level = 1-x)
#'    return( IC.tempo[["coef"]][3,"upper"])
#' }
#'
#' summary(fm1)$tTable
#' discreteRoot(fctIC2,grid = seq(1/1000,1,0.001), increasing = FALSE)$par
#' discreteRoot(fctIC3,grid = seq(1/1000,1,0.001), increasing = FALSE)$par
#'
#' ## negative coefficient
#' fctIC <- function(x){    
#'    IC.tempo <- intervals(fm1, level = x)
#'    return( IC.tempo[["coef"]][3,"upper"])
#' }
#' discreteRoot(fctIC,grid = seq(0,1-1/1000,0.001))$par

#' 
#' @export
discreteRoot <- function(fn, grid, increasing = TRUE, check = TRUE,
                         tol = .Machine$double.eps ^ 0.5) {

    n.grid <- length(grid)    
    value.grid <- rep(NA, n.grid)    
    iter <- 1
    ncv <- TRUE
    iSet <- 1:n.grid
    factor <- c(-1,1)[increasing+1]
    
### ** Check
    if(check){
        value.grid[1] <- fn(grid[1])
        value.grid[n.grid] <- fn(grid[n.grid])
        if(sign(value.grid[1])==value.grid[n.grid]){
            list(par = NA,
                 value = NA,
                 counts = 2,
                 cv = 1,
                 message = "Cannot find a solution because the function does not change sign \n")
        }
    }

    
### ** Expore the grid using dichotomic search
    while(iter <= n.grid && ncv && length(iSet)>0){
        iMiddle <- ceiling(length(iSet)/2)
        iIndexInSet <- iSet[iMiddle]
        if(check==FALSE || iIndexInSet %in% c(1,n.grid) == FALSE){
            value.grid[iIndexInSet] <- fn(grid[iIndexInSet])
        }        
        if(is.na(value.grid[iIndexInSet])){
            iSet <- setdiff(iSet,iMiddle)
            iter <- iter + 1
        }else if(factor*value.grid[iIndexInSet] > tol){
            iSet <- iSet[setdiff(1:iMiddle,iMiddle)]
            iter <- iter + 1
        }else if(factor*value.grid[iIndexInSet] < -tol){
            iN.set <- length(iSet)
            iSet <- iSet[setdiff(iMiddle:iN.set,iMiddle)]
            iter <- iter + 1
        }else{
            ncv <- FALSE
            solution <- grid[iIndexInSet]
            value <- value.grid[iIndexInSet]
        }
        
    }

### ** If did not find a value whose image matched tol, give the closest solution
    if(ncv){
       iIndexInSet <- which.min(abs(value.grid))

       ncv <- FALSE
       solution <- grid[iIndexInSet]
       value <- value.grid[iIndexInSet]
    }

    return(list(par = solution,
                value = value,
                counts = iter,
                cv = ncv,
                message = NULL))
}


##----------------------------------------------------------------------
### discreteRoot.R ends here
