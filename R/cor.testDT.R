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
  
  butils.base::validCharacter(format, validValues = c("wide","long"), validLength = 1, method = "cor.test.data.table")
  butils.base::validCharacter(output, validValues = c("data.table","matrix","plot"), validLength = 1, method = "cor.test.data.table")
  
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
    butils.base::validNames(names, validValues = names(data))
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

