#' @title Perform bootstrap computation on an object
#' 
#' @description Perform bootstrap computation on an object under H0 or H1. Handle one grouping variable.
#' 
#' @param object the fitted model.
#' @param data the data that have been used to fit the model.
#' @param n.boot the number of replications. Should be a large number.
#' @param n.cpus the number of cpu to use.
#' @param export.publish use \code{Publish::publish} to define the parameter of interest
#' @param FUN.coef the function used to extract the statistics of interest from the model.
#' @param var.id the variable in the dataset indicating the grouping level of the data. 
#' @param var.group the group variable use to resample under H0. Otherwise resampling is done under H1.
#' @param load.library additional library to load
#' @param seed set the random number generator
#' 
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' Bootstrap under H1: randomly select observations (or individuals according to argument var.id) to form a new dataset.
#' If the same individual appear several times, a different group value is given for each apparition.
#'
#' Bootstrap under H0: randomly permute the group label (according to argument var.group) at the observation level (or individuals according to argument var.id).
#' 
#' The Publish package must be installed for the function to work when setting the argument \code{export.publish} to \code{TRUE}.
#' The Publish package is available on Github: https://github.com/tagteam/Publish
#' and can be installed using, for e.g., the devtools package: devtools::install_github("tagteam/Publish")
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
#' resG <- bootGLS(m.lm, n.boot = 1e2, var.group = "G")
#' 
#' res
#' print(res, seq.length.out = 10)
#' resG 
#' }
#' 
#' 
#' \dontrun{ 
#' res <- bootGLS(m.lm, n.boot = 1e4, n.cpus = 4)
#' resG <- bootGLS(m.lm, n.boot = 1e4, var.group = "G", n.cpus = 4)
#' 
#' res
#' print(res, seq.length.out = 10)
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
#' bootGLS(fm1, n.boot = 1e2, var.id = "Mare")
#' }
#' \dontrun{
#' res <- bootGLS(fm1, n.boot = 1e4, n.cpus = 4, var.id = "Mare")
#' apply(res$all.boot, 2, mean)
#' }
#' 
#' }
#' 
#' @export
bootGLS <- function(object,
                    n.boot,
                    FUN.coef,
                    data = NULL,
                    n.cpus = 1,
                    export.publish = FALSE,
                    var.id = NULL,
                    var.group = NULL,
                    load.library = "nlme",
                    seed = 10){
  
  if(is.null(data)){
    data <- extractData(object, force = TRUE, convert2dt = TRUE)
  }else{
    data <- copy(as.data.table(data))
  }
    
  ## tests
  if(is.null(var.id)){
    if("idBOOT" %in% names(data)){
      stop("\"idBOOT\" must not be the name of any column in object$call$data \n")
    }
    data$idBOOT <- 1:NROW(data)
    var.id <- "idBOOT"
  }else if(var.id %in% names(data) == FALSE){
    stop("wrong specification of argument \'var.id\' \n",
         var.id," does not name any column in object$call$data \n",
         "possible names: ",paste(names(data), collapse = " "),"\n")
  }
  if(!is.null(var.group) && var.group %in% names(data) == FALSE){
    stop("wrong specification of argument \'var.group\' \n",
         var.group," does not name any column in object$call$data \n",
         "possible names: ",paste(names(data), collapse = " "),"\n")
  }
  
  ## publish
  if(export.publish){
    requireNamespace("Publish")
    object.publish <- Publish::publish(object, print = FALSE)  
    load.library <- c("Publish", load.library)
  }
  ## 
  setkeyv(data, var.id)
  vec.id <- unique(data[[var.id]]) #  all ids
  if(!is.null(var.group)){
    trueGroup <- data[,.SD[[1]][1], by = var.id, .SDcols = var.group][[2]]
  }else{
    indexId <- lapply(vec.id, function(x){which(data[[var.id]]==x)})
    n.obsId <- unlist(lapply(indexId, length))
  }
  
  ##
  if(export.publish){
    variable <- object.publish$rawTable[,"Variable"]
    n.coef <- length(variable)
    
    index.blanck <- which(variable == "")
    variable[index.blanck] <- variable[sapply(index.blanck, function(x){max(setdiff(1:x,index.blanck))})]
    e.coef <- setNames(object.publish$rawTable[,"Coefficient"],
                       paste(variable,object.publish$rawTable[,"Units"])
    )
  }else{
    e.coef <- coef(object)
    n.coef <- length(e.coef)
  }
  n.id <- length(vec.id)
  
  ## 
  if(export.publish){
    FUN.coef <- function(x){
      res <- Publish::publish(x, print = FALSE)
      return(res$rawTable$Coefficient)
    }
  }
  if(missing(FUN.coef)){
    if(class(object) == "lme"){
      FUN.coef <- "fixef"
    }else{
      FUN.coef <- "coef"
    }
  }
  
  ## boot
  if (!missing(seed)) set.seed(seed)
  bootseeds <- sample(1:max(1e6,seed),size=n.boot,replace=FALSE)
  
  resampleFCT <- function(){
    dt.new <- copy(data)
    index <- sample(n.id,replace=TRUE) # random of ids
    
    if(!is.null(var.group)){
      bootGroup <- trueGroup[index] # permutation of the group label attributing one of a random individual
      dt.new[, (var.group) := bootGroup[.GRP], by = var.id] # update of the group label at the individual level
    }else{
      dt.new <- dt.new[unlist(indexId[index])] # randomly pick individuals to form the new dataset
      newID <- unlist(lapply(1:length(index), function(i){ rep(i,n.obsId[index[i]])})) # form unique ID for individuals
      dt.new[, (var.id) := newID] # set the new ids
    }
    return(dt.new)
  }
  
  if(n.cpus>1){
    
    cl <- parallel::makeCluster(n.cpus)
    doParallel::registerDoParallel(cl)
    
    boots <- foreach::`%dopar%`(foreach::foreach(b=1:n.boot,.packages=c(load.library,"data.table"),.export=NULL), {
      set.seed(bootseeds[b])
      
      object$call$data <- resampleFCT()
      
      objectBoot <- tryCatch(do.call(FUN.coef,list(eval(object$call))),
                             error = function(x){return(NULL)})    
      return(objectBoot)
    })
    # perform bootstrap
    M.boot <- do.call(rbind, boots)
    
    if(is.null(M.boot)){ # failure of all bootstraps
      boots <- foreach::`%dopar%`(foreach::foreach(b=1:1,.packages=c(load.library,"data.table"),.export=NULL), {
        set.seed(bootseeds[b])
        object$call$data <- resampleFCT()
        do.call(FUN.coef,list(eval(object$call)))
      })
    }
    
  }else{
    
    M.boot <- matrix(NA, nrow = n.boot, ncol = n.coef)
    pb <- utils::txtProgressBar(max = n.boot)
    for(b in 1:n.boot){
      set.seed(bootseeds[b])
      
      object$call$data <- resampleFCT()
      objectBoot <- tryCatch(do.call(FUN.coef,list(eval(object$call))),
                             error = function(x){return(NULL)})    
      if(!is.null(objectBoot)){
        M.boot[b,] <- objectBoot
      }
      utils::setTxtProgressBar(pb, value = b)
      
    }
    close(pb)
    
  }
  
  ## post treatment
  CI <- calcCI(M.boot, var.group, e.coef, n.coef)
  p.value <- calcPvalue(M.boot, var.group, e.coef, n.coef)
  
  ## export
  ls.export <- list(all.boot = M.boot,
                    p.value = p.value,
                    coef = e.coef,
                    CI = CI,
                    var.group = var.group,
                    n.boot = n.boot)
  class(ls.export) <- "glsboot"
  return(ls.export)
}

# {{{ associated functions

#' @title Display the result fo the bootstrap computation
#' @description Display the result fo the bootstrap computation
#' 
#' @param x object obtained with the function \code{bootGLS}
#' @param seq.length.out the number of subsamples on which the results should be computed
#' @param ... not used
#'
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' @method print glsboot
#' @export
print.glsboot <- function(x, seq.length.out, ...){
  
  rowM <- rbind(estimate = x$coef,
                estimate.boot = apply(x$all.boot,2,median,na.rm = TRUE),
                x$CI,
                p.value = x$p.value)    
  print(t(rowM))
  n.bootReal <- colSums(!is.na(x$all.boot))
  cat("for ",as.double(min(n.bootReal))," replications of the resampling \n", sep = "")
  
  if(!missing(seq.length.out)){
    cat("\n")
    resSubset <- calcStatSubset(all.boot = x$all.boot, 
                                coef = x$coef, 
                                var.group = x$var.group, 
                                length.out = seq.length.out,
                                n.boot = x$n.boot)
    
    print(resSubset)
    cat("for subsets \n")
  }
  
  return(invisible(x))
}


calcStatSubset <- function(all.boot, coef, var.group, length.out, n.boot, n.coef){
  n.coef <- length(coef)
  
  seqIndex <- round(seq(1,n.boot+1, length.out = length.out+1))
  M.index <- cbind(seqIndex[-length(seqIndex)],(seqIndex-1)[-1])
  n.Index <- NROW(M.index)
  
  resSubset <- lapply(1:n.Index, function(i){
    M.boot_red <- all.boot[M.index[i,1]:M.index[i,2],,drop = FALSE]
    CI_tempo <- calcCI(M.boot_red, var.group, coef, n.coef)
    p.value_tempo <- calcPvalue(M.boot_red, var.group, coef, n.coef)
    return(rbind(CI_tempo,p.value = p.value_tempo))
  })
  
  names(resSubset) <- paste0(M.index[,1],":",M.index[,2])
  return(resSubset)
}

calcCI <- function(M.boot, var.group, e.coef, n.coef){
  
  if(!is.null(var.group)){
    names.coef <- names(e.coef)
    CI <- sapply(1:n.coef, function(row){
      if(names.coef[row]==var.group){
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

calcPvalue <- function(M.boot, var.group, e.coef, n.coef){
  
  if(!is.null(var.group)){
    names.coef <- names(e.coef)
    p.value <- sapply(1:n.coef, function(row){
      if(names.coef[row]==var.group){
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
  
  optimum <- optim(par = log(-log(0.5)), fn=function(p){
    p2 <- exp(-exp(p)) # convert to [0;1] scale
    abs(quantile(dist, probs = p2))
  }, method = "L-BFGS-B")
  optimum$par <- exp(-exp(optimum$par)) # backtransform
  
  if(optimum$par>0.5){
    p.value <- 2*(1-optimum$par)
  }else{
    p.value <- 2*optimum$par
  }
  
  return(p.value)
}

# }}}
