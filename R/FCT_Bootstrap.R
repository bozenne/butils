#' @title Perform bootstrap computation on an object
#' @param object the fitted model.
#' @param n.boot the number of replications. Should be a large number.
#' @param n.cpus the number of cpu to use.
#' @param IDvar the variable in the dataset indicating the grouping level of the data. 
#' @param GROUPvar the group variable use to resample under H0. Otherwise resampling is done under H1.
#' @param load.library additional library to load
#' 
#' @details
#' Should have an attribute call that includes a data attribute (i.e. object$call$data must exist).
#'
#' @examples
#' set.seed(10)
#' n <- 100
#' n.boot <- 1e4
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' df <- rbind(data.frame(Y=Y1,G=1),
#'            data.frame(Y=Y2,G=2)
#'            )
#' m.lm <- lm(Y ~ G, data = df)
#' bootGLS(m.lm, n.boot = 10, n.cpus = 2, GROUPvar = "G")

bootGLS <- function(object, n.boot, n.cpus, IDvar = NULL, GROUPvar = NULL,
                    load.library = "nlme", seed = 10){
  
  library(parallel)
  library(doParallel)
  library(data.table)
  
  data <- as.data.table(copy(eval(object$call$data)))
  
  ## tests
  if(is.null(IDvar)){
    if("idBOOT" %in% names(data)){
      stop("\"idBOOT\" must not be the name of any column in object$call$data \n")
    }
    data$idBOOT <- 1:NROW(data)
    IDvar <- "idBOOT"
  }else if(IDvar %in% names(data) == FALSE){
    stop("wrong specification of argument IDvar \n",
         IDvar," does not name any column in object$call$data \n",
         "possible names: ",paste(names(data), collapse = " "),"\n")
  }
  if(!is.null(GROUPvar) && GROUPvar %in% names(data) == FALSE){
    stop("wrong specification of argument GROUPvar \n",
         GROUPvar," does not name any column in object$call$data \n",
         "possible names: ",paste(names(data), collapse = " "),"\n")
  }
  
  ## 
  vec.id <- unique(data[[IDvar]]) #  all ids
  if(!is.null(GROUPvar)){
    trueGroup <- data[,.SD[[1]][1],by = IDvar, .SDcols = GROUPvar][[2]]
  }else{
    indexId <- lapply(vec.id, function(x){which(data[[IDvar]]==x)})
  }
  
  ##
  e.coef <- coef(object)
  n.coef <- length(e.coef)
  n.id <- length(vec.id)
  
  ## boot
  if (!missing(seed)) set.seed(seed)
  bootseeds <- sample(1:max(1e6,seed),size=seed,replace=FALSE)
  
  cl <- parallel::makeCluster(n.cpus)
  doParallel::registerDoParallel(cl)
  
  boots <- foreach(b=1:n.boot,.packages=c(load.library,"data.table"),.export=NULL) %dopar% {
    set.seed(bootseeds[[b]])
    
    # input dt.Neutral, vec.id, trueGroup, myModel
    dt.new <- copy(data)
    index <- sample(n.id,replace=FALSE) # random of ids
    
    if(!is.null(GROUPvar)){
      bootGroup <- trueGroup[index]
      dt.new[, group := bootGroup[.GRP], by = IDvar]
    }else{
      dt.new <- dt.new[unlist(indexId[index])]
    }
    
    object$call$data <- dt.new
    
    objectBoot <- tryCatch(eval(object$call),
                           error = function(x){return(NULL)})    
    
    return(coef(objectBoot))
  }
  
  #
  M.boot <- do.call(rbind, boots)
  if(is.null(GROUPvar)){
    
    p.value <- sapply(1:n.coef, function(row){
      mean( abs(M.boot[,row]) >= abs(e.coef[row]) )
    })
    CI <- sapply(1:n.coef, function(row){
      quantile(M.boot[,row], probs = c(0.025,0.975))+e.coef[row]
    })
    
  }else{
    
    p.value <- sapply(1:n.coef, function(row){
      findP1(M.boot[,row])
    })  
    CI <- sapply(1:n.coef, function(row){
      quantile(M.boot[,row], probs = c(0.025,0.975))
    })
    
  }
  colnames(CI) <- names(e.coef)
  names(p.value) <- names(e.coef)
  
  
  return(list(all.boot = M.boot,
              p.value = p.value,
              CI = CI,
              n.boot = n.boot))
}

findP1 <- function(dist){
  2*optim(par = 0.5, lower = 0, upper = 1, fn=function(p){
    abs(quantile(dist, probs = p))      
  }, method = "L-BFGS-B")$par
}
