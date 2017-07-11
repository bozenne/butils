### getSigmaGLS.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 26 2017 (13:56) 
## Version: 
## last-updated: maj 26 2017 (14:17) 
##           By: Brice Ozenne
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Variance-covariance of GLS object
#' 
#' @description Rebuilt the variance-covariance matrix from a GLS object
#' 
#' @param gls the gls object
#' @param type Can be either "correlation" or "covariance".
#' @param individual cluster for which the variance covariance matrix should be returned
#' @param upper logical value indicating whether the upper triangle of the distance matrix should be returned
#' @param addId the id name at the row and col names
#' @param plot should the correlation matrix be displayed. \emph{logical}.
#' @param args.plot a list of arguments to be passed to \code{ggHeatmap} to specify how the correlation matrix should be displayed
#' @param output how to output the correlation value. Can be \code{matrix}, \code{data.table} or \code{plot}.
#' @param trace should the progression of the computation of the correlation be displayed. \emph{logical}.
#' 
#' @examples 
#' if(require(nlme)){
#' data(Ovary)
#' Ovary <- as.data.table(Ovary)
#' setkeyv(Ovary, "Mare")
#' fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), 
#'             data = Ovary, correlation = corAR1(form = ~ 1 | Mare))
#' getSigmaGLS(fm1)
#' }
#' @keywords function correlation display gls
#' @export
getSigmaGLS <- function(gls, type = "covariance", upper = NULL, individual = NULL, addId = TRUE,
                        plot = TRUE, args.plot = list(), output = "matrix", trace = 1){
  
  butils.base::validCharacter(output, validValues = c("data.table","matrix","plot"), validLength = 1, method = "getSigmaGLS")
  butils.base::validCharacter(type, validValues = c("correlation","covariance"), validLength = 1, method = "getSigmaGLS")
  data <- nlme::getData(gls)
  
  #### rebuilt the matrix of correlation
  if (!is.null(gls$modelStruct$corStruct)) {
    
    ## check the ordering of the data
    corObj <- eval(gls$call$correlation)
    formulaCor <- get("formula.corStruct", asNamespace("nlme"))(corObj)
    variableCor <- all.vars(nlme::getGroupsFormula(formulaCor))
    
    if(any(order(data[[variableCor]]) != seq_len(NROW(data)))){
      stop("getSigmaGLS: the dataset used to fit the model must be sort by \"",paste(variableCor, collapse = "\" \""),"\"",
           " to ensure valid extraction of the variance covariance matrix \n")
    }
    if(is.numeric(data[[variableCor]])){
      stop("getSigmaGLS: the cluster variable \"",paste(variableCor, collapse = "\" \""),"\"",
           " must be a character/factor vector to ensure valid extraction of the variance covariance matrix \n")
    }
    
    ## find the cluster for which the varCov should be returned
    ls.Sigma <- nlme::corMatrix(gls$modelStruct$corStruct)
    Id <- names(ls.Sigma)
    if(is.null(individual)){
      if(!is.list(ls.Sigma)){
        individual <- Id[1]
      }else if(all(unlist(lapply(ls.Sigma, identical, ls.Sigma[[1]])) == TRUE)){
        individual <- Id[1]
      }else {
        seqRow <- unlist(lapply(ls.Sigma, nrow))
        individual <- Id[which.max(seqRow)]
        
        if(trace>0){
          cat("getSigmaGLS: variance-covariance matrices differ between clusters \n",
              "return the first largest matrix (\"",paste(variableCor, collapse = "\" \""),"\"=",individual,")\n",sep = "")
        }
      }
      
    }else if(is.numeric(individual)){
      individual <- Id[individual]
    }
    
    ## correlation matrix
    Sigma <- nlme::corMatrix(gls$modelStruct$corStruct)[individual]
    Sigma <- as.matrix(Matrix::bdiag(Sigma))
    
  } else {
    n.rep <- length(coef(gls$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE))
    Sigma <- diag(1,n.rep)
  }
  
  #### name matrix
  if (!is.null(individual) && !is.null(gls$modelStruct$corStruct)) {
    index.individual <- which(gls$groups %in% individual)
    nameCluster <- gls$groups[index.individual]
  }else{
    nameCluster <- NULL
  }
  
  if(!is.null(gls$modelStruct$varStruct)) {
    varObj <- eval(gls$call$weights)
    formulaVar <- get("formula.varFunc", asNamespace("nlme"))(varObj)
    variableVar <- paste(all.vars(nlme::getGroupsFormula(formulaVar)), collapse = " ")
    
    if(!is.null(individual) && !is.null(gls$modelStruct$corStruct)){
      nameRep <- names(nlme::varWeights(gls$modelStruct$varStruct)[index.individual])
    }else{
      nameRep <- names(coef(gls$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE))  
    }
  }else{
    names.time <- as.character(1:nrow(Sigma))
    variableVar <- "repetition"
    nameRep <- NULL
  }
  names.time <- paste(nameRep," (",nameCluster,")",sep="")
  colnames(Sigma) <- names.time
  rownames(Sigma) <- names.time
  
  #### add the variance terms
  if (type == "covariance"){
    if(!is.null(gls$modelStruct$varStruct)) {
      if(is.null(gls$groups)){
        vw <- 1/nlme::varWeights(gls$modelStruct$varStruct)[names.time]
      }else{
        ind <- gls$groups %in% individual ### need data to be sorted by individual
        vw <- 1/nlme::varWeights(gls$modelStruct$varStruct)[ind]
        names.time <- names(vw)
      }
      
    } else{ 
      vw <- rep(1, nrow(Sigma))
    }
    vars <- (gls$sigma * vw)^2
    Sigma <- t(Sigma * sqrt(vars)) * sqrt(vars)
  }
  
  ## manage matrix
  if(!is.null(upper)){
    if(upper){
      gdata::lowerTriangle(Sigma) <- NA
    }else{
      gdata::upperTriangle(Sigma) <- NA
    }
  }
  
  #### convert to data table
  index <- expand.grid(1:nrow(Sigma),1:nrow(Sigma))
  dt.Sigma <- data.table(covariance = as.numeric(t(Sigma)),
                         var1 = rownames(Sigma)[index[,1]],
                         var2 = colnames(Sigma)[index[,2]],
                         rep = cumsum(duplicated(rownames(Sigma)))
  )
  setnames(dt.Sigma, c("covariance","var1","var2"), c(type,paste0(variableVar,1),paste0(variableVar,2)))
 
  if(plot || output == "plot"){
    if("xlab" %in% names(args.plot) == FALSE){args.plot$xlab <- variableVar}
    if("ylab" %in% names(args.plot) == FALSE){args.plot$ylab <- variableVar}
    if("legend_title" %in% names(args.plot) == FALSE){args.plot$legend_title <- type}
    
  
    dt.Sigma2 <- copy(dt.Sigma)
    dt.Sigma2[[paste0(variableVar,1)]] <- make.unique(rownames(Sigma))[index[,1]]
    dt.Sigma2[[paste0(variableVar,2)]] <- make.unique(rownames(Sigma))[index[,2]]
    
    gg <- do.call("ggHeatmap", args = c(list(data = dt.Sigma2), 
                                        list(name.x = paste0(variableVar,1), 
                                             name.y = paste0(variableVar,2),
                                             name.fill = type), 
                                        args.plot))
    if(output != "plot"){print(gg)}
  }
  
  #### export
  if(output == "data.table"){
    return(dt.Sigma)
  }else if(output == "matrix"){
    return(Sigma)
  }else if(output == "plot"){
    return(gg)
  }
  
}


#----------------------------------------------------------------------
### getSigmaGLS.R ends here
