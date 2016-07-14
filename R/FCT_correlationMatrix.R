#' @title Perform pairwise correlation tests between matrix columns
#' 
#' @export
cor.testDT <- function(data, names = NULL, format, lower.tri = TRUE,
                       method.cor = "cor.test", use.pairwiseNNA = TRUE, 
                       reorder = "AOE", imput.value = 0, hclus.method = "complete",
                       trace = TRUE, plot = TRUE, args.plot = list(), output = "data.table", ...){
  
  validCharacter(format, validValues = c("wide","long"), validLength = 1, method = "cor.test.data.table")
  validCharacter(output, validValues = c("data.table","matrix","plot"), validLength = 1, method = "cor.test.data.table")
  
  ## data format
  if(is.data.table(data)){
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
    dataL <- melt(data, id.vars = NULL, measure.vars = names, value.name = "value", variable.name = "variable")
    name.variable <- "variable"
    name.value <- "value"
  }else{
    name.variable <- names[1]
    name.value <- names[2]
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
    order <- corrplot::corrMatOrder(Wcorr, order = reorder, hclust.method = hclus.method)
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
                                   dt <- as.data.table(reshape2::melt(cor.array[,,x], value.name = dimnames(cor.array)[[3]][x]))
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
  if(output == "data.table"){
    return(cor.dt)
  }else if(output == "matrix"){
    return(cor.array)
  }else if(output == "plot"){
    return(gg)
  }
  
}

#' @export
ggHeatmap <- function(data, name.x, name.y, name.fill, add.text, round = NULL,
                      title = "", xlab = "", ylab = "",  legend_title = "correlation",
                      na.value = "grey50", col_low = "blue", col_midpoint = "white", col_high = "red", midpoint = 0, limits = NULL, 
                      textSize = 15, angle.x = 90){
  
  if(is.data.table(data)){
    data <- as.data.table(data)
  }else{
    data <- copy(data)
  }
  
  if(is.null(limits)){
    limits <- range(data[[name.fill]], na.rm = TRUE)
  }
  
  gg <- ggplot(data, aes_string(x = name.x, y = name.y, fill = name.fill)) + geom_tile()
  gg <- gg + ggtitle(title)
  gg <- gg + xlab(xlab) + ylab(ylab)
  gg <- gg + scale_fill_gradient2(low = col_low, mid = col_midpoint, high = col_high, 
                                  name = legend_title, 
                                  midpoint = midpoint, na.value = na.value, limits = limits)
  gg <- gg + theme(text = element_text(size=textSize), axis.text.x = element_text(angle = angle.x, hjust = 1))
  
  if(!missing(add.text)){
    
    if(!is.null(round)){
      if(is.numeric(round)){
        data[, add.text := round(.SD[[1]], digit = round), .SDcols = add.text]
      }else if(round == "p.value"){
        type <- findInterval(data[["p.value"]], vec = c(0.001,0.01,0.05,0.1))
        type <- factor(type, levels = 0:4, labels = c("***","**","*",".",""))
        data[, p.value := type]
      }else{
        stop("ggHeatmap: non valid value for argument \'round\' \n")
      }
    }
    
    gg <- gg +  geom_text(aes_string(fill = name.fill, label = add.text))  
  }
  
  return(gg)
}