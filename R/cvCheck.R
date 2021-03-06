## * cvCheck
#' @title Test the sensibility of the lvm estimate to the initialization points
#' @description Test the sensibility of the lvm estimate to the initialization points
#' @name cvCheck
#'
#' @param object a lvm model
#' @param data a data frame
#' @param factor.vcov inflation factor for the variance when sampling the initialization points
#' @param n.init number of initialization points to be used
#' @param ncpus the number of CPU to be used
#' @param keep.cov should the covariance between parameter be kept to simulate the initialization points
#' @param trace should a progression bar be displayed?
#' @param ... additional arguments to be passed to estimate
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples 
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' covariance(m) <- v1~v2+v3+v4
#' dd <- lava::sim(m,10000) ## Simulate 10000 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' \dontrun{
#' summary(cvCheck(m, dd, ncpus = 1))
#' # summary(cvCheck(m, dd, ncpus = 4))
#' }
#' \dontshow{
#' summary(cvCheck(m, dd, ncpus = 1, n.init = 3))
#' }
#'
#' 
#' @export
cvCheck <- function (object, ...) {
  UseMethod("cvCheck", object)
}

## * cvCheck.lvm
#' @rdname cvCheck
#' @export
cvCheck.lvm <- function(object,
                        data,
                        factor.vcov = 1,
                        n.init = 100,
                        keep.cov = TRUE,
                        ncpus = 1,
                        trace = TRUE,
                        ...){
    
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}

    dots <- list(...)
    if("control" %in% names(dots) == FALSE){
        dots$control <- list()
    }
   
    ## ** automatic intialisation
    test.W <- 0
    
    if(trace){cat("* initialisation \n")}
    suppressWarnings(
        lvm.init <- estimate(object, data, ...)  
    )

    e.coef <- stats::coef(lvm.init)
    name.coef <- names(e.coef)
    n.coef <- length(name.coef)
    var.param <- grep(",",name.coef, fixed = TRUE)
    mean.param <- setdiff(1:n.coef, var.param)

    lvm.vcov <- stats::vcov(lvm.init)
    if(keep.cov){        
        VCOV <- lvm.vcov
    }else{
        VCOV <- matrix(0, nrow = NROW(lvm.vcov), ncol = NCOL(lvm.vcov))        
    }
    diag(VCOV) <-  diag(lvm.vcov)*factor.vcov

    requireNamespace("tmvtnorm")
    sample.start <- tmvtnorm::rtmvnorm(n = n.init, mean = e.coef, sigma = VCOV,
                                       lower = c(rep(-Inf, length(mean.param)),
                                                 rep(0, length = length(var.param))),
                                       algorithm = "gibbs"
                                       )
    
    ## ** warper
    warper <- function(x){
        dots$control$start <- sample.start[x,]

        suppressWarnings(
            resStart <- try(do.call("estimate", 
                                    c(list(x = object, data = data), dots))
                            )
        )
        
        if(class(resStart) != "try-error"){
            return(c(resStart$opt$convergence,
                     as.numeric(stats::logLik(resStart)),
                     stats::coef(resStart) ))
        }else{
            return(rep(NA, n.coef+2))
        }
    }

    ## ** parallel computations
    if(ncpus>1){
      cl <- parallel::makeCluster(ncpus)
      doSNOW::registerDoSNOW(cl)
      
      if(trace > 0){
        pb <- utils::txtProgressBar(max = n.init, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
      }else{
        opts <- NULL
      }
      
      i <- NULL # [:for CRAN check] (foreach)
      Mres <- foreach::`%dopar%`(
        foreach::foreach(i = 1:n.init,
                         .packages =  "lava",
                         .export = "dots",
                         .combine = "cbind",
                         .options.snow = opts),
        {                                
          return(warper(i))
        })
      
      if(trace > 0){  close(pb) }
      parallel::stopCluster(cl)
    }else{
      if(trace>0){
        requireNamespace("pbapply")
        resLoop <- pbapply::pblapply(1:n.init, warper)
      }else{
        resLoop <- lapply(1:n.init, warper)
      }
      Mres <- do.call("cbind",resLoop)
    }

    ## ** postprocess and export
    Mres <- cbind(Mres,
                  c(lvm.init$opt$convergence,
                    as.numeric(stats::logLik(lvm.init)),
                    stats::coef(lvm.init) ))  
    df.resCV <- stats::setNames(as.data.frame(t(Mres)), c("cv", "logLik",name.coef))
    
    ## ** export
    out <- list(estimates = df.resCV,
                seeds = sample.start)
    class(out) <- c("cvlvm","data.frame")
    return(out)
}

## * summary.cvlvm
#' @title Summary function associated with cvCheck
#' @description Summary function associated with cvCheck
#'
#' @param object the output of cvCheck
#' @param threshold threshold used to cluster the convergence points (height in hclust)
#' @param ... for compatibility only
#'
#' @method summary cvlvm
#' @export
summary.cvlvm <- function(object, threshold = NULL, ...){


    e <- object$estimates
    index.cv <- which(e$cv==0)    
    n.init <- length(e$cv)
    
    quantile.cv <- apply(e[index.cv,-1], 2, function(x){
        c(range = diff(range(x, na.rm = TRUE)), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
    })
    cat("Summary of the estimates: \n", sep = "")
    print(quantile.cv)
    

    pc.cv <- mean(e$cv==0, na.rm = TRUE)
    cat("Convergence rate: ",round(100*pc.cv,2),"% (",sum(e$cv==0, na.rm = TRUE)," over ",n.init,") \n", sep = "")

    if(!is.null(threshold)>0 && pc.cv>0){
    n.cv <- length(index.cv)
    indexRed.bestcv <- which.max(e$logLik[index.cv])
            
    dist.coef <- stats::dist(e[index.cv,-(1:2),drop=FALSE]) 
    hclust.res <- stats::hclust(dist.coef, method="complete")
    n.clusters <- 1 + sum(hclust.res$height > threshold)
    
    dist.coef <- sum(as.matrix(dist.coef)[indexRed.bestcv,] < threshold)
    
    cat("Threshold: ",threshold," \n", sep = "")
    cat("Number of convergence points: ",n.clusters," \n", sep = "")
    cat("Convergence rate to the best convergence point: ",round(100*dist.coef/n.cv,2),"% (",dist.coef," over ",n.cv,") \n", sep = "")
    print(e[index.cv[indexRed.bestcv],])
    
  }else{
    quantile.cv <- NA
    n.clusters <- NA
  }
  
  return(invisible(list(quantile = quantile.cv,
                        pc.cv = pc.cv,
                        n.clusters = n.clusters)))
}


## * optimx1
optimx1 <- function(start,objective,gradient,hessian,...) {
  requireNamespace("optimx")
    
  nlminbcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min","optim.method","optimx.method")
  dots <- list(...)
  control <- list(...)$control
  control <- control[names(control)%in%nlminbcontrols]
  dots$control <- control
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  mypar <- c(list(start=start,objective=objective,gradient=gradient,hessian=hessian),dots)
  mypar["debug"] <- NULL
  
  if(!is.null(mypar$control$optimx.method)){
    optimx.method <- mypar$control$optimx.method
    mypar$control$optimx.method <- NULL
    
    optimx.res <- optimx::optimx(par = mypar$start, 
                                 fn = mypar$objective, 
                                 gr = mypar$gradient,
                                 method = optimx.method,
                                 control = mypar$control)
    
  }else{
    mypar$control$all.methods <- TRUE
    
    optimx.res <- optimx::optimx(par = mypar$start, 
                                 fn = mypar$objective, 
                                 gr = mypar$gradient,
                                 control = mypar$control)  
  }
  
  
  index.cv <- which(optimx.res$convcode == 0)
  if(length(index.cv) == 0){
    res <- list(par = rep(NA, length(mypar$start)),
                objective = NA,
                convergence = 1,
                evaluations = c("function" = NA, "gradient" = NA),
                message = "optimx has not converged \n"
    )
  }else if(length(index.cv) == 1){
    
    par_tempo <- unlist(lapply(1:length(mypar$start), function(x){optimx.res[[x]][index.cv]})) 
    names(par_tempo) <- names(mypar$start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cv],
                convergence = optimx.res$convcode[index.cv],
                evaluations = c("function" = optimx.res$fevals[index.cv], "gradient" = optimx.res$gevals[index.cv]),
                message = paste0("One algorithm has converged: ",attributes(optimx.res)$row.names[index.cv], " (value = ",optimx.res$value[index.cv],")\n")
    )
  }else{
    index.cvBest <- index.cv[which.max(optimx.res$value[index.cv])]
    par_tempo <- unlist(lapply(1:length(mypar$start), function(x){optimx.res[[x]][index.cvBest]})) 
    names(par_tempo) <- names(mypar$start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cvBest],
                convergence = optimx.res$convcode[index.cvBest],
                evaluations = c("function" = optimx.res$fevals[index.cvBest], "gradient" = optimx.res$gevals[index.cvBest]),
                message = paste0("Several algorithms have converged but only the best convergence point is retained \n",
                                 "method: ",paste(attributes(optimx.res)$row.names[index.cv], collapse = " | "), "\n",
                                 "value : ",paste(optimx.res$value[index.cv], collapse = " | "), "\n")
                
    )
  }
  
  return(res)
  
}
