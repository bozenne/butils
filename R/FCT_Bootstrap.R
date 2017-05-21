#' @title Perform bootstrap computation on an object
#' 
#' @description Perform bootstrap computation on an object under H0 or H1. Handle one grouping variable.
#' 
#' @param object the fitted model.
#' @param n.boot the number of replications. Should be a large number.
#' @param n.cpus the number of cpu to use.
#' @param fctCoef the function used to extract the statistics of interest from the model.
#' @param IDvar the variable in the dataset indicating the grouping level of the data. 
#' @param GROUPvar the group variable use to resample under H0. Otherwise resampling is done under H1.
#' @param load.library additional library to load
#' @param seed set the random number generator
#' 
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' Bootstrap under H1: randomly select observations (or individuals according to argument IDvar) to form a new dataset.
#' If the same individual appear several times, a different group value is given for each apparition.
#'
#' Bootstrap under H0: randomly permute the group label (according to argument GROUPvar) at the observation level (or individuals according to argument IDvar).
#' 
#' @examples
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
#' res <- bootGLS(m.lm, n.boot = 1e4, n.cpus = 4)
#' resG <- bootGLS(m.lm, n.boot = 1e4, GROUPvar = "G", n.cpus = 4)
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
#' bootGLS(fm1, n.boot = 1e2, IDvar = "Mare")
#' }
#' \dontrun{
#' res <- bootGLS(fm1, n.boot = 1e4, n.cpus = 4, IDvar = "Mare")
#' apply(res$all.boot, 2, mean)
#' }
#' 
#' }
#' 
#' @export
bootGLS <- function(object,
                    n.boot,
                    fctCoef,
                    n.cpus = 1,
                    IDvar = NULL,
                    GROUPvar = NULL,
                    load.library = "nlme",
                    seed = 10){
  data <- try(as.data.table(copy(eval(object$call$data))), silent = TRUE)
  if("try-error" %in% class(data)){
    data <- as.data.table(copy(eval(object$call$data, envir = parent.frame())))
  }
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
    setkeyv(data, IDvar)
    vec.id <- unique(data[[IDvar]]) #  all ids
    if(!is.null(GROUPvar)){
        trueGroup <- data[,.SD[[1]][1], by = IDvar, .SDcols = GROUPvar][[2]]
    }else{
        indexId <- lapply(vec.id, function(x){which(data[[IDvar]]==x)})
        n.obsId <- unlist(lapply(indexId, length))
    }
    
    ##
    e.coef <- coef(object)
    n.coef <- length(e.coef)
    n.id <- length(vec.id)
      
    ## 
    if(missing(fctCoef)){
      if(class(object) == "lme"){
        fctCoef <- "fixef"
      }else{
        fctCoef <- "coef"
      }
    }
    
    ## boot
    if (!missing(seed)) set.seed(seed)
    bootseeds <- sample(1:max(1e6,seed),size=n.boot,replace=FALSE)

    resampleFCT <- function(){
      dt.new <- copy(data)
      index <- sample(n.id,replace=TRUE) # random of ids
      
      if(!is.null(GROUPvar)){
        bootGroup <- trueGroup[index] # permutation of the group label attributing one of a random individual
        dt.new[, (GROUPvar) := bootGroup[.GRP], by = IDvar] # update of the group label at the individual level
      }else{
        dt.new <- dt.new[unlist(indexId[index])] # randomly pick individuals to form the new dataset
        newID <- unlist(lapply(1:length(index), function(i){ rep(i,n.obsId[index[i]])})) # form unique ID for individuals
        dt.new[, (IDvar) := newID] # set the new ids
      }
      return(dt.new)
    }
    
    if(n.cpus>1){
      
      cl <- parallel::makeCluster(n.cpus)
      doParallel::registerDoParallel(cl)
      
      boots <- foreach::`%dopar%`(foreach::foreach(b=1:n.boot,.packages=c(load.library,"data.table"),.export=NULL), {
        set.seed(bootseeds[b])
        
        object$call$data <- resampleFCT()
        
        objectBoot <- tryCatch(do.call(fctCoef,list(eval(object$call))),
                               error = function(x){return(NULL)})    
        
        return(objectBoot)
      })
      
      # perform bootstrap
      M.boot <- do.call(rbind, boots)
      
    }else{
      
      M.boot <- matrix(NA, nrow = n.boot, ncol = n.coef)
      pb <- utils::txtProgressBar(max = n.boot)
      for(b in 1:n.boot){
        set.seed(bootseeds[b])
        
        object$call$data <- resampleFCT()
        
        objectBoot <- tryCatch(do.call(fctCoef,list(eval(object$call))),
                               error = function(x){return(NULL)})    
        if(!is.null(objectBoot)){
          M.boot[b,] <- objectBoot
        }
        utils::setTxtProgressBar(pb, value = b)
        
      }
      close(pb)
      
    }
    
    ## post treatment
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
#' @description Display the result fo the bootstrap computation
#' 
#' @param x object obtained with the function \code{bootGLS}
#' @param seq_length.out the number of subsamples on which the results should be computed
#' @param ... not used
#'
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' @method print glsboot
#' @export
print.glsboot <- function(x, seq_length.out, ...){

    rowM <- rbind(estimate = x$coef,
                  estimate.boot = apply(x$all.boot,2,median,na.rm = TRUE),
                  x$CI,
                  p.value = x$p.value)    
    print(t(rowM))
    n.bootReal <- colSums(!is.na(x$all.boot))
    cat("for ",as.double(min(n.bootReal))," replications of the resampling \n", sep = "")

    if(!missing(seq_length.out)){
        cat("\n")
      resSubset <- calcStatSubset(all.boot = x$all.boot, 
                     coef = x$coef, 
                     GROUPvar = x$GROUPvar, 
                     length.out = seq_length.out,
                     n.boot = x$n.boot)
       
        print(resSubset)
        cat("for subsets \n")
     }

    return(invisible(x))
}


calcStatSubset <- function(all.boot, coef, GROUPvar, length.out, n.boot, n.coef){
  n.coef <- length(coef)
  
  seqIndex <- round(seq(1,n.boot+1, length.out = length.out+1))
  M.index <- cbind(seqIndex[-length(seqIndex)],(seqIndex-1)[-1])
  n.Index <- NROW(M.index)
  
  resSubset <- lapply(1:n.Index, function(i){
    M.boot_red <- all.boot[M.index[i,1]:M.index[i,2],,drop = FALSE]
    CI_tempo <- calcCI(M.boot_red, GROUPvar, coef, n.coef)
    p.value_tempo <- calcPvalue(M.boot_red, GROUPvar, coef, n.coef)
    return(rbind(CI_tempo,p.value = p.value_tempo))
  })
  
  names(resSubset) <- paste0(M.index[,1],":",M.index[,2])
  return(resSubset)
}

calcCI <- function(M.boot, GROUPvar, e.coef, n.coef){

      if(!is.null(GROUPvar)){
        names.coef <- names(e.coef)
        CI <- sapply(1:n.coef, function(row){
            if(names.coef[row]==GROUPvar){
                quantile(M.boot[,row], probs = c(0.025,0.975), na.rm = TRUE)+e.coef[row]
            }else{
                return(c("2.5%"=NA,"97.5%"=NA))
            }
        })

    }else{

        CI <- sapply(1:n.coef, function(row){
            quantile(M.boot[,row], probs = c(0.025,0.975), na.rm = TRUE)
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
                mean( abs(M.boot[,row]) >= abs(e.coef[row]), na.rm = TRUE)
            }else{
                return(NA)
            }
        })
        
    }else{
        p.value <- sapply(1:n.coef, function(row){
            findP1(na.omit(M.boot[,row]))
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
