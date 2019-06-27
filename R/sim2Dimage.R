### sim2Dimage.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2018 (09:39) 
## Version: 
## Last-Updated: jun 27 2019 (09:40) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * sim2Dimage (documentation)
##' @title Simulate a 2D Image
##' @description Simulate a 2D image with a given mean and autocorrelation structure.
##' @name sim2Dimage
##'
##' @param n [integer, >0] number of images to be simulated
##' @param coords [2 column data.frame] coordinates of the element in the image.
##' @param x.sites [integer, >0] number of x-coordinates. Ignored when \code{coords} is specified.
##' @param y.sites [integer, >0] number of y-coordinates. Ignored when \code{coords} is specified.
##' @param mu [numeric vector] Expected value at each site.
##' @param vgm [variogramModel] Variogram used to generate the background noise. Output of the \code{vgm} function.
##' @param name.X [character] name of the column for the design matrix.
##'
##' @return A list with two elements \itemize{
##' \item coords: the coordinates of the sites
##' \item X: the design matrix containing the signal of each site for each replicate.
##' }


## * sim2Dimage (examples)
##' @rdname sim2Dimage
##' @examples
##'
##' if(require(sp) & require(gstat)){
##' out0 <- sim2Dimage(3, x.sites = 50, y.sites = 50,
##'                   mu = 0,
##'                  vgm = vgm(psill = 1, range = 2, model='Exp'))
##' dim(out0$X)
##' spdf0 <- cbind(out0$coords,out0$X[1,])
##' gridded(spdf0) <- ~x+y
##' plot(spdf0)
##'
##' coords <- expand.grid(x = 1:50,
##'                       y = 1:50)
##' coords$mu <- 0
##' coords$mu[coords$x<=5 | coords$y<=5] <- 2
##'
##' out1 <- sim2Dimage(3, coords = coords[,c("x","y")],
##'                   mu = coords$mu,
##'                  vgm = vgm(psill = 1, range = 2, model='Exp'))
##' spdf1 <- cbind(out1$coords,out1$X[1,])
##' gridded(spdf1) <- ~x+y
##' plot(spdf1)
##'
##' }

## * sim2Dimage (code)
##' @rdname sim2Dimage
##' @export
sim2Dimage <- function(n, coords = NULL, x.sites, y.sites, mu, vgm, name.X = "X"){


    ## ** spatial coordinates
    if(is.null(coords)){
        coords <- expand.grid(x = 1:x.sites,
                              y = 1:y.sites)
        coords.init <- TRUE
    }else{
        coords.init <- FALSE
    }
    n.sites <- NROW(coords)
    if(length(mu) == 1){
        mu <- rep(mu, n.sites)
    }
    if(length(mu) != n.sites){
        if(coords.init){
            stop("Length of argument \'mu\' must be equal to \'x.sites\' times \'y.sites\' \n")
        }else{
            stop("Length of argument \'mu\' must be equal to the length of argument \'coords\' \n")
        }
    }
    
    ## ** simulate background noise and add mean
    requireNamespace("gstat")
    g.model <- gstat::gstat(formula = z~1,
                            locations = ~x+y,
                            dummy = T,
                            beta = 0,
                            model = vgm,
                            nmax = 20)
    X <- matrix(NA, nrow = n, ncol = n.sites)
    for(iN in 1:n){ ## iN <- 1
        mess <- utils::capture.output(res <- predict(g.model, newdata = coords, nsim = 1, trace = FALSE))
        X[iN,] <- res[,3] + mu 
    }
    colnames(X) <- paste0(name.X,1:n.sites)
    
    return(list(coords = coords,
                X = X))

}


######################################################################
### sim2Dimage.R ends here
