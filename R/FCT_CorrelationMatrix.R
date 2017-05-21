#' @title Perform pairwise correlation tests between matrix columns
#' 
#' @description Compute pairwise correlation between multiple outcomes. Can reorder the outcomes according to a clustering of the obtained correlations.
#' 
#' @param data the dataset. 
#' @param format the format of data. Can be \code{"wide"} or \code{"long"} 
#' @param names the names of the column in the dataset to be analysed. If \code{NULL}, all the columns will be used.
#' @param lower.tri should only the lower traingle of the matrix be filled?
#' @param method.cor the method used to compute the correlation. 
#' @param use.pairwiseNNA If \code{FALSE} correlation is set to NA if their is one NA among the two outcomes. Otherwise the correlation is computed on pairs without NA.
#' @param reorder argument order of corrMatOrder. If \code{NULL} the original order of the variables is kept.
#' @param imput.value the value to be used in the clustering algorithm in place of NA. Can be a number or an operator (e.g. median). 
#' @param hclust.method argument hclust.method of corrMatOrder.
#' @param trace should the progression of the computation of the correlation be displayed. \emph{logical}.
#' @param plot should the correlation matrix be displayed. \emph{logical}.
#' @param args.plot a list of arguments to be passed to \code{ggHeatmap} to specify how the correlation matrix should be displayed
#' @param output how to output the correlation value. Can be \code{matrix}, \code{data.table} or \code{plot}.
#' @param ... additional arguments to be passed to method.cor
#' 
#' @details 
#' data must be coercible to data.table. \cr
#' \cr
#' If \code{data="wide"} each column correspond to a different outcomes, indicated by argument names. \cr
#' If \code{data="long"} then the first value of names is used to define the outcomes and the seconde is used to define their values. \cr
#' \cr
#' method.cor can be any method but its first two arguments, named x and y, will be given the values of the outcomes. 
#' \cr
#' If imput.value is an operator then it is applied on the values in the same line or column.
#' \cr 
#' This function uses the function \code{ggHeatmap} to display the correlation matrix.
#' 
#' @return An object determined by the output argument
#' 
#' @seealso \code{\link{ggHeatmap}} to display the correlation matrix
#' 
#' @examples 
#' M <- matrix(rnorm(1e3),100,10)
#' cor.testDT(M, format = "wide")
#' 
#' @keywords function correlation test
#' 
#' @export
cor.testDT <- function(data, format, names = NULL, lower.tri = TRUE,
                       method.cor = "cor.test", use.pairwiseNNA = TRUE, 
                       reorder = "AOE", imput.value = 0, hclust.method = "complete",
                       trace = TRUE, plot = TRUE, args.plot = list(), output = "data.table", ...){
  
  validCharacter(format, validValues = c("wide","long"), validLength = 1, method = "cor.test.data.table")
  validCharacter(output, validValues = c("data.table","matrix","plot"), validLength = 1, method = "cor.test.data.table")
  
  ## data format
  if(is.data.table(data) == FALSE){
    data <- as.data.table(data)
  }else{
    data <- copy(data)
  }
  
  ## reduce data
  if(is.null(names)){
    names <- names(data)
  }else{
    validNames(names, validValues = names(data))
  }
  n.name <- length(names)
  
  ## reshape data
  n <- nrow(data)
  if(format == "wide"){
    dataL <- data.table::melt(data, id.vars = NULL, measure.vars = names, value.name = "value", variable.name = "variable")
    name.variable <- "variable"
    name.value <- "value"
  }else{
    name.variable <- names[1]
    name.value <- names[2]
    dataL <- data
  }
  setkeyv(dataL, name.variable)
  
  ## compute the correlation
  cor.array <- array(NA, dim = c(n.name, n.name, 6), dimnames = list(names, names, c("n.NNA", "n.NA", "correlation", "CIinf", "CIsup", "p.value")))
  if(trace){ pb <- utils::txtProgressBar(max = n.name*n.name) }
  
  for(iterName1 in 1:n.name){
    for(iterName2 in iterName1:n.name){
      col1 <- dataL[names[iterName1]][[name.value]]
      col2 <- dataL[names[iterName2]][[name.value]]
      index.NNA <- intersect(which(!is.na(col1)), which(!is.na(col2)))
      n.NNA <- length(index.NNA)
      
      if((n == n.NNA) || (n.NNA>0 && use.pairwiseNNA) ){
        test <- do.call(method.cor, args = list(x = col1[index.NNA], y = col2[index.NNA], ...)) 
        
        newValues <- list(n.NNA = n.NNA, n.NA = n - n.NNA,
                          correlation = test$estimate, CIinf = test$conf.int[1], CIsup = test$conf.int[2], 
                          p.value = test$p.value)
        if(is.null(newValues$CIinf)){newValues$CIinf <- NA}
        if(is.null(newValues$CIsup)){newValues$CIsup <- NA}
        
      }else{
        newValues <- c(n.NNA = n.NNA, n.NA = n - n.NNA, 
                       correlation = NA, CIinf = NA, CIsup = NA, 
                       p.value = NA)
      }
      # cor.dt <- rbind(cor.dt, c(variable1 =  names[iterName2], variable2 = names[iterName1], newValues))
      cor.array[names[iterName2],names[iterName1],] <- unlist(newValues)
      cor.array[names[iterName1],names[iterName2],] <- unlist(newValues)
      # cor.dt <- rbind(cor.dt, c(variable1 =  names[iterName1], variable2 = names[iterName2], newValues))
      
      if(trace){ utils::setTxtProgressBar(pb = pb, value = (iterName1-1)*n.name+iterName2)}
    } 
  }
  
  if(trace){ close(pb) }
  
  #### reordering
  Wcorr <- cor.array[,,"correlation"]
  if(!is.null(reorder)){
    
    if(any(is.na(Wcorr))){ #### imputation for missing values
      if(is.numeric(imput.value)){
        Wcorr[is.na(Wcorr)]  <- imput.value
      }else{
        index.NA <- which(is.na(Wcorr), arr.ind = TRUE)
        value.NA <- apply(index.NA, 1, function(x){
          do.call(imput.value, list(na.omit(c(Wcorr[x[1],],Wcorr[,x[2]]))))
        })
        Wcorr[index.NA] <- value.NA
      }
    }  
    order <- corrplot::corrMatOrder(Wcorr, order = reorder, hclust.method = hclust.method)
    cor.array <- cor.array[order,order,]
  }
  
  #### remove useless values
  if(lower.tri == FALSE){
    for(iter3 in 1:dim(cor.array)[3]){
      cor.array[,,iter3][upper.tri(cor.array[,,iter3])] <- NA
    }
  }
  
  #### convert to data.table
  cor.dt <- Reduce(cbind, lapply(1:6, 
                                 function(x){
                                   dt <- as.data.table(data.table::melt(cor.array[,,x], value.name = dimnames(cor.array)[[3]][x]))
                                   if(x>1){dt[, c("Var1","Var2") := NULL]}
                                   return(dt)
                                 }))
  setnames(cor.dt, old = c("Var1", "Var2"), new = paste0(name.variable,1:2))
  
  #### display
  if(plot || output == "plot"){
    gg <- do.call("ggHeatmap", args = c(list(data = cor.dt), 
                                        list(name.x = paste0(name.variable,1)), 
                                        list(name.y = paste0(name.variable,2)), 
                                        list(name.fill = "correlation"), 
                                        args.plot))
    if(output != "plot"){print(gg)}
  }
  
  #### output
  out <- list()
  if(output == "data.table"){
    out$dt <- cor.dt
  }
  if(output == "matrix"){
    out$array <- cor.array
  }
  if(output == "plot"){
    out$plot <- gg
  }
  
  return(out)
}

#' @title Display pairwise correlation
#' 
#' @description Display a correlation matrix (or heatmap) using ggplot2
#' 
#' @param data the dataset containing the correlation values.
#' @param name.x name of the column containing the label of the first outcome involved in the correlation
#' @param name.y name of the column containing the label of the second outcome involved in the correlation
#' @param name.fill name of the column containing the value of the correlation
#' @param add.text additional information to be displayed above the squares representing the correlation. Must name a column in data.
#' @param round if not \code{NULL}, the number of digit used to round add.text using signif. Can also be \code{p.value} to indicates with the usual convention the significance level.
#' @param title the title of the plot
#' @param xlab the name of the x axis
#' @param ylab the name of the y axis
#' @param legend_title the name of the legend
#' @param na.value the color used to display missing values (NA).
#' @param col_low the color used to diplay low correlation values
#' @param col_midpoint the color used to diplay intermediate correlation values
#' @param col_high the color used to diplay high correlation values
#' @param midpoint the value corresponding to intermediate correlation value.
#' @param limits the minimum and maximum correlation values used to build the color panel
#' @param textSize the size of the text in the plot
#' @param angle.x the inclination of the x labels
#' 
#' @details 
#' data must be coercible to data.table. \cr
#' 
#' @return a ggplot object
#' 
#' @keywords function correlation display
#' @export
ggHeatmap <- function(data, name.x, name.y, name.fill, add.text, round = NULL,
                      title = "", xlab = "", ylab = "",  legend_title = "correlation",
                      na.value = "grey50", col_low = "blue", col_midpoint = "white", col_high = "red", midpoint = 0, limits = NULL, 
                      textSize = 15, angle.x = 90){
  
  if(!is.data.table(data)){
    data <- as.data.table(data)
  }else{
    data <- copy(data)
  }
  
  #### prepare
  if(!missing(add.text)){
    
    if(add.text %in% names(data) == FALSE){
      stop("ggHeatmap: variable ",add.text," not found \n",
           "add.text should be one of \"",paste(names(data), collapse = "\" \""),"\"\n")
    }
    
    if(!is.null(round)){
      if(is.numeric(round)){
        data[, add.text := round(.SD[[1]], digit = round), .SDcols = add.text]
      }else if(round == "p.value"){
        type <- findInterval(data[["p.value"]], vec = c(0.001,0.01,0.05,0.1))
        type <- factor(type, levels = 0:4, labels = c("***","**","*",".",""))
        data[, add.text := type]
      }else{
        stop("ggHeatmap: non valid value for argument \'round\' \n")
      }
    }else{
      data[, add.text := .SD[[1]], .SDcols = add.text]
    }
    
  }
  
  data[[name.x]] <- factor(data[[name.x]], levels = unique(data[[name.x]]))
  data[[name.y]] <- factor(data[[name.y]], levels = unique(data[[name.y]]))
  
  if(is.null(limits)){
    limits <- range(data[[name.fill]], na.rm = TRUE)
  }
  
  #### plot
  gg <- ggplot(data, aes_string(x = name.x, y = name.y, fill = name.fill)) + geom_tile()
  gg <- gg + ggtitle(title)
  gg <- gg + xlab(xlab) + ylab(ylab) +  scale_y_discrete(limits = rev(levels(data[[name.x]])))
  gg <- gg + scale_fill_gradient2(low = col_low, mid = col_midpoint, high = col_high, 
                                  name = legend_title, 
                                  midpoint = midpoint, na.value = na.value, limits = limits)
  gg <- gg + theme(text = element_text(size=textSize), axis.text.x = element_text(angle = angle.x, hjust = 1))
  
  if(!missing(add.text)){
    gg <- gg +  geom_text(aes_string(fill = name.fill, label = "add.text"))  
  }
  
  return(gg)
}


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
#' Ovary <- Ovary[order(Ovary$Mare),]
#' fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), 
#'             data = Ovary, correlation = corAR1(form = ~ 1 | Mare))
#' getSigmaGLS(fm1)
#' }
#' @keywords function correlation display gls
#' @export
getSigmaGLS <- function(gls, type = "covariance", upper = NULL, individual = NULL, addId = TRUE,
                        plot = TRUE, args.plot = list(), output = "matrix", trace = 1){
  
  validCharacter(output, validValues = c("data.table","matrix","plot"), validLength = 1, method = "getSigmaGLS")
  validCharacter(type, validValues = c("correlation","covariance"), validLength = 1, method = "getSigmaGLS")
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
