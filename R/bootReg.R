## * bootReg - documentation
#' @title Perform bootstrap computation on an object
#' @name bootReg
#' 
#' @description Perform bootstrap computation on an object under H0 or H1. Handle one grouping variable.
#' 
#' @param object the fitted model.
#' @param data the data that have been used to fit the model.
#' @param strata if not \code{NULL}, a stratified bootstrap is performed according to this variable.
#' @param type the type of test for which the bootstrap should be performed. Can be \code{"coef"}, \code{"anova"}, \code{"publish"}.
#' Setting type to \code{NULL} enable the use of \code{FUN.estimate} and \code{FUN.stdError}.
#' @param cluster the variable indicating the level where the sample is i.i.d. Only required for gls with no correlation argument.
#' @param n.boot the number of replications. Should be a large number.
#' @param n.cpus the number of cpu to use.
#' @param FUN.estimate the function used to extract the punctual estimates from the model.
#' @param FUN.stdError the function used to extract the standard error associated with the punctual estimate (i.e. standard error of the empirical estimator).
#' @param FUN.resample the function used simulate new data under the model. Default is \code{NULL} which corresponds to a non-parametric bootstrap.
#' @param FUN.iid the function used to extract the influence function from the model.
#' @param load.library additional library to load on each CPU. Useful when performing parallel computation.
#' @param seed set the random number generator
#' @param trace should the execution of the bootstrap be displayed using a progress bar?
#' @param name.cluster internal argument.
#' @param rejectIfWarning Should the estimate be ignored if a warning is returned by the estimation routine?
#' @param ... ignored
#' 
#' @details
#' Bootstrap under H1: randomly select observations (or individuals according to argument var.id) to form a new dataset.
#' If the same individual appear several times, a different group value is given for each apparition.
#' 
#' @examples
#' #### data  ####
#' n <- 1e2
#' set.seed(10)
#' df.data <- data.frame(Y = rnorm(n),
#'                      group = gl(3, 5, n, labels = c("Ctl","Trt","Neu")),
#'                      gender = gl(2, 5, n, labels = c("Female","Male"))[sample.int(n)]
#'                      )
#' 
#' #### lm ####
#' m.lm <- lm(Y ~ group*gender, data = df.data)
#' resBoot <- bootReg(m.lm, n.boot = 1e2)
#' resBoot
#' summary(resBoot, type = "norm")
#' summary(resBoot, type = "basic")
#' summary(resBoot, type = "stud")
#' summary(resBoot, type = "perc")
#' summary(resBoot, type = "bca")
#' 
#' resBoot <- bootReg(m.lm, FUN.resample = "simulate", n.boot = 1e1)
#' resBoot
#' 
#' #### gls ####
#' library(nlme)
#' e.gls <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
#'              data = Ovary, correlation = corAR1(form = ~ 1 | Mare))
#' resBoot <- bootReg(e.gls, n.boot = 1e1)
#'
#' #### lme ####
#' e.lme <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
#'              data = Ovary, random =~ 1 | Mare)
#' resBoot <- bootReg(e.lme, n.boot = 1e1)
#' 
#' @export
`bootReg` <-
  function(object, ...) UseMethod("bootReg")


## * bootReg.lm
#' @rdname bootReg
#' @export
bootReg.lm <- function(object,
                       type = "coef",
                       FUN.estimate = NULL,
                       FUN.stdError = NULL,
                       data = NULL,
                       load.library = NULL,
                       ...){

### ** extract data
    if(is.null(data)){
        data <- as.data.table(extractData(object, model.frame = FALSE))
    }else{
        data <- copy(as.data.table(data))
    }

### ** add cluster variable
    if("XXclusterXX" %in% names(data)){
        stop("\"XXclusterXX\" must not be the name of any column in data \n")
    }
    data$XXclusterXX <- 1:NROW(data)
    name.cluster <- "XXclusterXX"

### ** define default statistic

    if(!is.null(FUN.estimate)){
        
    }else if(identical(type,"coef")){
        FUN.stdError <- function(x){
            return(summary(x)$coef[,"Std. Error"])
        }        
        FUN.estimate <- function(x){
            return(summary(x)$coef[,"Estimate"])
        }
    }else if(identical(type,"anova")){
        FUN.estimate <- function(x){ # x <- object
            res <- anova(x)            
            return(setNames(res[,"F value"],rownames(res)))
        }        
    }else if(identical(type,"publish")){
        requireNamespace("Publish")
        load.library <- c("Publish", load.library)        
        FUN.estimate <- function(x){ # x <- object
            res <- Publish::publish(x, print = FALSE)
            return(setNames(res$rawTable$Coefficient, res$rawTable$Variable))
        }
        
    }else{
        stop("arguments \'FUN.estimate\' must be specified \n",
             "when type is not \"coef\" \"anova\" \"publish\" \n")
    }

### ** run boostratp
    out <- .bootReg(object,
                    data = data, name.cluster = name.cluster,
                    FUN.estimate = FUN.estimate,
                    FUN.stdError = FUN.stdError,
                    load.library = load.library,
                    ...)

### ** export
    return(out)
    
}

## * bootReg.gls
#' @rdname bootReg
#' @export
bootReg.gls <- function(object,
                        cluster,
                        type = "coef",
                        FUN.estimate = NULL,
                        FUN.stdError = NULL,
                        data = NULL,                        
                        load.library = "nlme",
                        ...){

### ** extract data
    if(is.null(data)){
        data <- as.data.table(extractData(object, model.frame = FALSE))
    }else{
        data <- copy(as.data.table(data))
    }

### ** identify cluster variable
    cluster2 <- as.numeric( nlme::getGroups(object))
    
    if(length(cluster2)!=0){ ## no correlation
        if(NCOL(cluster2)==1){            
            cluster <- colnames(object$groups)
            if(is.null(cluster)){
                cluster <- all.vars(nlme::getGroupsFormula(object$modelStruct$corStruct))
            }
        }else{
            stop("Can only handle one grouping variable \n")
        }
    }else{
        test.invalid <- length(cluster) != 1 || !is.character(cluster) || (cluster %in% names(data) == FALSE)
        if(test.invalid){
            stop("Argument \'cluster\' must refer to a column in data when there is no groupping variable defined by the model \n")
        }
    }

    ### ** define default statistic

    if(!is.null(FUN.estimate)){
        
    }else if(identical(type,"coef")){
        FUN.estimate <- function(x){
            return(summary(x)$tTable[,"Value"])
        }
        FUN.stdError <- function(x){
            return(summary(x)$tTable[,"Std.Error"])
        }
    }else if(identical(type,"anova")){
        FUN.estimate <- function(x){ # x <- object
            res <- anova(x)            
            return(setNames(res[,"F-value"],rownames(res)))
        }      
    }else if(identical(type,"publish")){
        requireNamespace("Publish")
        load.library <- c("Publish", load.library)
        
        FUN.estimate <- function(x){ # x <- object
            res <- Publish::publish(x, print = FALSE)
            return(setNames(res$rawTable$Coefficient, res$rawTable$Variable))
        }        
    }else{
        stop("arguments \'FUN.estimate\' must be specified \n",
             "when type is not \"coef\" \"anova\" \"publish\" \n")
    }

### ** run boostratp
    out <- .bootReg(object,
                    data = data, name.cluster = cluster,
                    FUN.estimate = FUN.estimate,
                    FUN.stdError = FUN.stdError,
                    load.library = load.library,
                    ...)

### ** export
    return(out)
    
}

## * bootReg.lme
#' @rdname bootReg
#' @export
bootReg.lme <- bootReg.gls

## * bootReg.lvmfit
#' @rdname bootReg
#' @export
bootReg.lvmfit <- function(object,
                           type = "coef",
                           FUN.estimate = NULL,
                           FUN.stdError = NULL,
                           data = NULL,
                           load.library = "lava",
                           ...){    
    
    ### ** extract data
    if(is.null(data)){
        data <- as.data.table(extractData(object, model.frame = FALSE))
    }else{
        data <- copy(as.data.table(data))
    }

    ### ** add cluster variable
    if("XXclusterXX" %in% names(data)){
        stop("\"XXclusterXX\" must not be the name of any column in data \n")
    }
    data$XXclusterXX <- 1:NROW(data)
    name.cluster <- "XXclusterXX"

    ### ** define default statistic
    if(!is.null(FUN.estimate)){
        
    }else if(identical(type,"coef")){
        FUN.stdError <- function(x){
            return(summary(x)$coef[names(x$opt$estimate),"Std. Error"])
        }        
        FUN.estimate <- function(x){
            return(summary(x)$coef[names(x$opt$estimate),"Estimate"])
        }
    }else{
        stop("arguments \'FUN.estimate\' must be specified \n",
             "when type is not \"coef\" \n")
    }

    ### ** run boostratp    
    out <- .bootReg(object,
                    data = data, name.cluster = name.cluster,
                    FUN.estimate = FUN.estimate,
                    FUN.stdError = FUN.stdError,
                    load.library = load.library,
                    ...)

### ** export
    return(out)
    
}

## * .bootReg function
#' @rdname bootReg
#' @export
.bootReg <- function(object, data, strata = NULL, name.cluster,
                     FUN.estimate, FUN.stdError, FUN.resample = NULL, FUN.iid = NULL,
                     n.boot = 1e3, n.cpus = 1, load.library, seed = 1, rejectIfWarning = TRUE,

                     trace = TRUE){

    .Random.seed_save <- .Random.seed
    
### ** extract coefficients
    e.coef <- do.call(FUN.estimate, args = list(object))
    if(!is.null(FUN.stdError)){
        e.stdError <- do.call(FUN.stdError, args = list(object))
    }else{
        e.stdError <- NULL
    }
    n.coef <- length(e.coef)
    name.coef <- names(e.coef)
    
### ** identify clusters
    setkeyv(data, name.cluster)
    unique.cluster = unique(data[[name.cluster]])
    n.cluster <- length(unique.cluster)
    
    ls.index.cluster <- lapply(unique.cluster, function(x){
        which(data[[name.cluster]]==x)
    })
    names(ls.index.cluster) <- unique.cluster
    n.obs.cluster <- unlist(lapply(ls.index.cluster,length))

    
    ### ** seed
    if (!is.null(seed)){
        set.seed(seed)
        bootseeds <- sample(1:max(1e6,seed),size=n.boot,replace=FALSE)
    }else{
        bootseeds <- NULL
    }
    
    ### ** non-parametric bootstrap simulation
    FUN.resample.save <- FUN.resample
    response.var <- NULL
    if(!is.null(strata)){
        if(any(strata %in% names(data) == FALSE)){
            stop("argument \'strata\' contains names that are not in data \n")
        }
    }
    
    if(is.null(FUN.resample)){
        FUN.resample <- function(object, data, response.var, ...){
            if(is.null(strata)){
                new.cluster <- sample(n.cluster,replace=TRUE) # random of ids
            }else{
                ## unique(.SD[[2]])   : contains the different ids of the patients (e.g. 1:10)
                ## .SD[[1]][1] (or .N): contains the number of rows (e.g. 10)
                ## sample(x = .N, replace=TRUE): randomly draw integers (e.g. 2, 5, 10, ...)
                res.sample <- data[,list(list(new.cluster = unique(.SD[[name.cluster]])[sample(x = .N, replace=TRUE)])),
                                   by=strata, .SDcols = name.cluster]
                new.cluster <- unlist(res.sample[[2]])
            }
            
            data.new <- data[unlist(ls.index.cluster[new.cluster])] # randomly pick individuals to form the new dataset
            newID <- unlist(lapply(1:n.cluster, function(i){ # i <- 1
                rep(i,n.obs.cluster[new.cluster[i]])
            })) # form unique ID for individuals
            data.new[, (name.cluster) := newID] # set the new ids
            return(data.new)
        }
    }else if(identical(FUN.resample,"simulate")){
        response.var <- all.vars(formula(object))[1]
        
        FUN.resample <- function(object, data, response.var, ...){
            data.new <- copy(data)
            data.new[[response.var]] <- as.double(unlist(stats::simulate(object)))
            return(data.new)
        }
    }

    ### ** boostrap function
    FUN.boostrap <- function(iB){
        if(!is.null(bootseeds)){
            set.seed(bootseeds[iB])
        }
        
        object$call$data <- do.call(FUN.resample, args = list(object = object, data = data, response.var = response.var))

        objectBoot <- matrix(NA, nrow = 2, ncol = n.coef,
                             dimnames = list(c("estimate","stdError"),name.coef)
                             )
        res <- tryWithWarnings(eval(object$call))
        if(!is.null(res$value)){            
            estimate.tempo <- do.call(FUN.estimate, args = list(res$value))
            objectBoot["estimate",names(estimate.tempo)] <- estimate.tempo
            if(!is.null(FUN.stdError)){
                stdError.tempo <- do.call(FUN.stdError, args = list(res$value))
                objectBoot["stdError",names(stdError.tempo)] <- stdError.tempo
            }
        }
        if(!is.null(res$warnings)){
            attr(objectBoot,"warning") <- res$warnings
        }
        if(!is.null(res$error)){
            attr(objectBoot,"error") <- res$error
        }
        return(objectBoot)
    }
  
### ** boot
  if(n.cpus>1){
      vec.export <- c("object","tryWithWarnings")
      b <- NULL ## for CRAN check
      
      if(trace){
          cl <- snow::makeSOCKcluster(n.cpus)
          doSNOW::registerDoSNOW(cl)
          pb <- txtProgressBar(min=1, max=n.boot, style=3)
          progress <- function(n){ setTxtProgressBar(pb, n) }
          
          ls.boot <- foreach::`%dopar%`(foreach::foreach(b=1:n.boot,.packages=c(load.library,"data.table"),
                                                         .options.snow = list(progress = progress), .export = vec.export), {          
                                                             return(FUN.boostrap(b))
                                                         })
          
      }else{
          cl <- parallel::makeCluster(n.cpus)
          doParallel::registerDoParallel(cl)
          ls.boot <- foreach::`%dopar%`(foreach::foreach(b=1:n.boot,.packages=c(load.library,"data.table"),
                                                         .export = vec.export), {          
              return(FUN.boostrap(b))
          })
      }

      stopCluster(cl)
      
      
  }else{

      if(trace){
          requireNamespace("pbapply")
          ls.boot <- pbapply::pblapply(1:n.boot, function(iB){
              return(FUN.boostrap(iB))
          })      
      }else{
          ls.boot <- lapply(1:n.boot, function(iB){
              return(FUN.boostrap(iB))
          })      
      }  
      
  }

    ### ** postprocess results
    index.ok <- which(unlist(lapply(ls.boot, function(iB){
        is.null(attr(iB,"error"))
    })))
    if(rejectIfWarning){
        index.ok <- intersect(index.ok,
                              which(unlist(lapply(ls.boot, function(iB){
                                  is.null(attr(iB,"warning"))
                              })))
                              )
    }

    Mestimate.boot <- do.call("rbind",lapply(ls.boot[index.ok], function(iB){iB["estimate",]}))
    if(!is.null(FUN.stdError)){
        MstdError.boot <- do.call("rbind",lapply(ls.boot[index.ok], function(iB){iB["stdError",]}))
    }else{
        MstdError.boot <- NULL
    }

    if(!is.null(FUN.iid)){
        iid <- do.call(FUN.iid, args = list(object))
    }else{
        iid <- NULL
    }
    
### ** export
    if(is.null(strata)){
        strata <- rep(1,NROW(data))
    }else{
        strata <- data[,interaction(.SD), .SDcols = strata]
    }

    if(!is.null(seed)){
        assign(x = ".Random.seed",
               value = .Random.seed_save,
               envir = globalenv())
    }
    
    out <- list(call = object$call,
                estimate = e.coef,
                stdError = e.stdError,
                iid = iid,
                boot.estimate = Mestimate.boot,
                boot.stdError = MstdError.boot,
                n.boot = n.boot,
                data = data,
                seed = seed,
                strata = strata,
                .Random.seed = .Random.seed_save,
                FUN.estimate = FUN.estimate,
                FUN.stdError = FUN.stdError,
                FUN.resample = FUN.resample.save,
                FUN.iid = FUN.iid,
                name.cluster = name.cluster
                )

    class(out) <- "bootReg"
    return(out)
}


