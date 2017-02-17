#' @title Perform bootstrap computation on an object
#' 
#' @description Perform bootstrap computation on an object
#' 
#' @param object the fitted model.
#' @param n.boot the number of replications. Should be a large number.
#' @param n.cpus the number of cpu to use.
#' @param IDvar the variable in the dataset indicating the grouping level of the data. 
#' @param GROUPvar the group variable use to resample under H0. Otherwise resampling is done under H1.
#' @param load.library additional library to load
#' @param seed set the random number generator
#' 
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' @examples
#' 
#' #### lm ####
#' set.seed(10)
#' n <- 100
#' n.boot <- 1e4
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' df <- rbind(data.frame(Y=Y1,G=1),
#'            data.frame(Y=Y2,G=2)
#'            )
#' m.lm <- lm(Y ~ G, data = df)
#' summary(m.lm)
#' 
#' \dontshow{
#' res <- bootGLS(m.lm, n.boot = 1e2)
#' resG <- bootGLS(m.lm, n.boot = 1e2, GROUPvar = "G")
#' 
#' res
#' print(res, seq_length.out = 10)
#' resG 
#' }
#' 
#' 
#' \dontrun{ 
#' res <- bootGLS(m.lm, n.boot = 1e4)
#' resG <- bootGLS(m.lm, n.boot = 1e4, GROUPvar = "G")
#' 
#' res
#' print(res, seq_length.out = 10)
#' resG
#' }
#' 
#' #### gls ####
#' if(require("nlme")){
#' data(Ovary)
#' fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
#'              correlation = corAR1(form = ~ 1 | Mare))
#' summary(fm1)
#' 
#' \dontshow{
#' bootGLS(fm1, n.boot = 1e2, n.cpus = 2, IDvar = "Mare")
#' }
#' \dontrun{
#' bootGLS(fm1, n.boot = 1e4, n.cpus = 2, IDvar = "Mare")
#' }
#' 
#' }
#' 
#' @export
bootGLS <- function(object,
                    n.boot,
                    n.cpus = 1,
                    IDvar = NULL,
                    GROUPvar = NULL,
                    load.library = "nlme",
                    seed = 10){
    
    data <- as.data.table(copy(eval(object$call$data)))
    
    ## tests
    if(is.null(IDvar)){
        if("idBOOT" %in% names(data)){
            stop("\"idBOOT\" must not be the name of any column in object$call$data \n")
        }
        data$idBOOT <- 1:NROW(data)
        IDvar <- "idBOOT"
    }else if(IDvar %in% names(data) == FALSE){
        stop("wrong specification of argument \'IDvar\' \n",
             IDvar," does not name any column in object$call$data \n",
             "possible names: ",paste(names(data), collapse = " "),"\n")
    }
    if(!is.null(GROUPvar) && GROUPvar %in% names(data) == FALSE){
        stop("wrong specification of argument \'GROUPvar\' \n",
             GROUPvar," does not name any column in object$call$data \n",
             "possible names: ",paste(names(data), collapse = " "),"\n")
    }
    
    ## 
    vec.id <- unique(data[[IDvar]]) #  all ids
    if(!is.null(GROUPvar)){
        trueGroup <- data[,.SD[[1]][1], by = IDvar, .SDcols = GROUPvar][[2]]
    }else{
        indexId <- lapply(vec.id, function(x){which(data[[IDvar]]==x)})
    }
    
    ##
    e.coef <- coef(object)
    n.coef <- length(e.coef)
    n.id <- length(vec.id)
    
    ## boot
    if (!missing(seed)) set.seed(seed)
    bootseeds <- sample(1:max(1e6,seed),size=n.boot,replace=FALSE)
    
    cl <- parallel::makeCluster(n.cpus)
    doParallel::registerDoParallel(cl)

    boots <- foreach::`%dopar%`(foreach::foreach(b=1:n.boot,.packages=c(load.library,"data.table"),.export=NULL), {
        set.seed(bootseeds[b])
        
        # input dt.Neutral, vec.id, trueGroup, myModel
        dt.new <- copy(data)
        index <- sample(n.id,replace=TRUE) # random of ids
        
        if(!is.null(GROUPvar)){
            bootGroup <- trueGroup[index]
            dt.new[, (GROUPvar) := bootGroup[.GRP], by = IDvar]
        }else{
            dt.new <- dt.new[unlist(indexId[index])]
        }
        
        object$call$data <- dt.new
        
        objectBoot <- tryCatch(eval(object$call),
                               error = function(x){return(NULL)})    

        # coef(objectBoot)
        return(coef(objectBoot))
    })
    
    #
    M.boot <- do.call(rbind, boots)

    CI <- calcCI(M.boot, GROUPvar, e.coef, n.coef)
    p.value <- calcPvalue(M.boot, GROUPvar, e.coef, n.coef)

    ## export
    ls.export <- list(all.boot = M.boot,
                      p.value = p.value,
                      coef = e.coef,
                      CI = CI,
                      GROUPvar = GROUPvar,
                      n.boot = n.boot)
    class(ls.export) <- "glsboot"
    return(ls.export)
}

# {{{ associated functions

#' @title Display the result fo the bootstrap computation
#' 
#' @description Display the result fo the bootstrap computation
#' 
#' @param x object obtained with the function \code{bootGLS}
#' @param seq_length.out the number of subsamples on which the results should be computed
#' 
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
print.glsboot <- function(x, seq_length.out){

    rowM <- rbind(estimate = x$e.coef,
                  x$CI,
                  p.value = x$p.value)    
    print(t(rowM))
    cat("for ",x$n.boot," replications of the resampling \n", sep = "")

    if(!missing(seq_length.out)){
        cat("\n")
       seqIndex <- seq(1,x$n.boot, length.out = seq_length.out)
        n.Index <- length(seqIndex)
        n.coef <- length(x$p.value)
        resSubset <- lapply(1:(n.Index-1), function(i){
            M.boot_red <- x$all.boot[seqIndex[i]:seqIndex[i+1],,drop = FALSE]
            CI_tempo <- calcCI(M.boot_red, x$GROUPvar, x$coef, n.coef)
            p.value_tempo <- calcPvalue(M.boot_red, x$GROUPvar, x$coef, n.coef)
            return(rbind(CI_tempo,p.value = p.value_tempo))
        })

        names(resSubset) <- paste0(seqIndex[-n.Index],":",seqIndex[-1])
        print(resSubset)
        cat("for subsets \n")
     }

    return(invisible(x))
}


calcCI <- function(M.boot, GROUPvar, e.coef, n.coef){

    if(!is.null(GROUPvar)){
        names.coef <- names(e.coef)
        CI <- sapply(1:n.coef, function(row){
            if(names.coef[row]==GROUPvar){
                quantile(M.boot[,row], probs = c(0.025,0.975))+e.coef[row]
            }else{
                return(c("2.5%"=NA,"97.5%"=NA))
            }
        })

    }else{

        CI <- sapply(1:n.coef, function(row){
            quantile(M.boot[,row], probs = c(0.025,0.975))
        })
    }

    colnames(CI) <- names(e.coef)
    return(CI)
}

calcPvalue <- function(M.boot, GROUPvar, e.coef, n.coef){

    if(!is.null(GROUPvar)){
        names.coef <- names(e.coef)
        p.value <- sapply(1:n.coef, function(row){
            if(names.coef[row]==GROUPvar){
                mean( abs(M.boot[,row]) >= abs(e.coef[row]) )
            }else{
                return(NA)
            }
        })
        
    }else{
        p.value <- sapply(1:n.coef, function(row){
            findP1(M.boot[,row])
        })  
        
    }

    names(p.value) <- names(e.coef)
    return(p.value)

}

findP1 <- function(dist){

    if(all(dist>0) || all(dist<0)){return(0)}
    
  optimum <- optim(par = 0.5, lower = 0, upper = 1, fn=function(p){
    abs(quantile(dist, probs = p))
  }, method = "L-BFGS-B")
  
  if(optimum$par>0.5){
      p.value <- 2*(1-optimum$par)
  }else{
      p.value <- 2*optimum$par
  }

  return(p.value)
}

# }}}
