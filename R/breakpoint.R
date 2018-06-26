### breakpoint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 26 2018 (09:13) 
## Version: 
## Last-Updated: jun 26 2018 (09:23) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * breakpoint (documentation)
#' @title Find One or Two Breakpoints
#' @description Find one or two breakpoints.
#' @name breakpoint
#'
#' @param object a \code{lm} object under the assumption of a linear relationship
#' @param pattern [vector of character] the number and type of breakpoints to be search. 0 indicates a flat line. 
#' @param var.breakpoint [character] the variable regarding which the breakpoints should be found.
#' @param init.breakpoint [list] a list containing initialization values for the breakpoints for each pattern.
#' @param n.iter [integer, >0] the maximum number of iterations used to estimates the breakpoints.
#' @param tol [numeric, >0] the maximum accpetable difference between two consecutive estimates of the breakpoints.
#' When reached, the estimation algorithm stops.
#' @param n.points [integer, >0] the number of points used to display the fit.
#' @param trace [0,1,2] trace the execution of the function.
#' @param digits [integer] how to round values that are displayed in the terminal.
#'
#'
#' @details
#' Argument pattern: 111 corresponding to three lines with different slopes
#' while 101 corresponds to a three lines where the middle one has null slope.
#' 
#' @references Muggeo, V. M. R. Estimating regression models with unknown break-points.
#' Statistics in medicine 2003; 22:3055-3071. 

## * breakpoint (example)
#' @rdname breakpoint
#' @examples
#' library(lava)
#' library(data.table)
#' library(ggplot2)
#'
#' #### simulate data
#' m <- lvm(Y1[0:0.1] ~ X1,
#'          Y2[1:0.1] ~ 0,
#'          X2 ~ 1,
#'          Y3[3:0.1] ~ -1*X3)
#'
#' distribution(m,~X1) <- uniform.lvm(0, 1)
#' distribution(m,~X2) <- uniform.lvm(1, 2)
#' distribution(m,~X3) <- uniform.lvm(2, 3)
#'
#' set.seed(10)
#' dt <- as.data.table(lava::sim(m, n = 1e2))
#' dtL <- melt(dt,
#'      measure.vars = list(paste0("X",1:3), paste0("Y",1:3)),
#'      value.name = c("X","Y"))
#' 
#' gg <- ggplot(dtL, aes(x = X))
#' gg <- gg + geom_point(aes(y = Y, color = variable))
#' gg
#'
#'
#' #### fit breakpoint regression
#' e.lm <- lm(Y~X, data = dtL)
#' resBP <- breakpoint(e.lm)
#' 
#' gg + geom_point(data = resBP$BP101$fit, aes(y = fit))
#'
#'
#' #### example from the package segmented
#' if(require(segmented)){
#' GS <- segmented(e.lm, psi = c(0.5,1.5))
#' cbind(value = resBP$BP111$value,
#'       se = resBP$BP111$se)
#' GS$psi
#' }
#' 
#' if(require(gridExtra)){
#' grid.arrange(resBP[["BP111"]]$plot + theme(text = element_text(size = 10), legend.position="none"),
#'              resBP[["BP101"]]$plot + theme(text = element_text(size = 10), legend.position="none"),
#'              resBP[["BP11"]]$plot + theme(text = element_text(size = 10), legend.position="none"),
#'              resBP[["BP10"]]$plot + theme(text = element_text(size = 10), legend.position="none"),
#'              resBP[["BP01"]]$plot + theme(text = element_text(size = 10), legend.position="none"), nrow = 1)
#' }
#'  

## * breakpoint (code)
#' @rdname breakpoint
#' @export
breakpoint <- function(object, pattern = c("111","101","10","11","01"), var.breakpoint = NULL, init.breakpoint = NULL,
                       n.iter = 10, tol = 1e-3, n.points = 1e3,
                       trace = 2, digits = -log10(tol)){


    ## ** find data
    dt <- model.frame(object)
    if(!is.data.table(dt)){
        data.table::setDT(dt)
    }else{
        dt <- data.table::copy(dt)
    }

    ## ** check arguments
    ## *** data
    reserved.names <- c("fit","Us","Vs","Us2","Vs2","gamma","beta")
    if(any(names(dt) %in% reserved.names)){
        txt <- names(dt[names(dt) %in% reserved.names])
        stop("data contains reserved names: \"",paste0(txt, collapse = "\" \""),"\"\n")
    }


    
    ## *** var.breakpoint
    if(is.null(var.breakpoint)){
        var.breakpoint <- all.vars(update(formula(object),"0~."))
    }

    if(length(var.breakpoint)!=1){
        stop("Argument \"var.breakpoint\" should correspond to only one variable \n")
    }

    if(var.breakpoint %in% names(dt) == FALSE){
        stop("Variable \"",var.breakpoint,"\" is not in data \n")
    }

    ## *** pattern
    if(any(pattern %in% c("111","101","10","11","01") == FALSE) ){
        stop("n.breakpoint must be one of \"111\",\"101\",\"10\",\"11\",\"01\" \n")
    }

    ## *** init.breakpoint
    if(is.null(init.breakpoint)){
        init.breakpoint <- list()
    }else{
        if(!is.null(init.breakpoint)){
            stop("Argument \'init.breakpoint\' must be a list \n")
        }
        if(any(pattern %in% names(init.breakpoint) == FALSE)){
            stop("Argument \'init.breakpoint\' must be a list containing elements \"",paste0(pattern,collpase = "\" \""),"\" \n")
        }
    }

    ## ** compute breakpoints
    for(iPattern in pattern){
        if(trace>0){ cat("* Search breaking points (pattern ",iPattern,")", sep = "") }
        if(trace>1){ cat("\n") }
        
        object[[paste0("BP",iPattern)]] <- .BPfit(object = object, data = dt, var.breakpoint = var.breakpoint, init.breakpoint = init.breakpoint[[iPattern]], pattern = iPattern,
                                                  n.iter = n.iter, tol = tol, n.points = n.points,
                                                  trace = (trace-1)>0, digits = digits)
        if(trace == 0){
            cat(" - done \n")
        }
        ## object$BP111$plot
    }
    
    ## ** export
    return(object)
}

## * .calcBP
.BPfit <- function(object, data, var.breakpoint, init.breakpoint, pattern,
                   n.iter, tol, n.points,
                   trace, digits){

    n.breakpoint <- nchar(pattern)-1
    
    if(is.null(init.breakpoint)){
        probs.breakpoint <- seq(0,1, length.out = n.breakpoint+2)[2:(n.breakpoint+1)]
        init.breakpoint <- quantile(data[[var.breakpoint]], probs = probs.breakpoint)
    }else if(length(init.breakpoint) != n.breakpoint){
        stop("Incorrect initialization for the breakpoints (too many or to few) \n")
    }

    ## ** lvm model
    ## free flat free
    var.response <- all.vars(update(formula(object),".~0"))
    if(pattern == "111"){
        formula.updated <- update(formula(object),".~.+Us+Vs+Us2+Vs2")
        
        coef.alpha <- var.breakpoint
        coef.beta <- c("Us","Us2")
        coef.gamma <- c("Vs","Vs2")
    }else if(pattern == "101"){
        txt.formula <- paste0(var.response,"~alpha*",var.breakpoint,"+beta*Us+Vs+Us2+Vs2")
        m <- lava::lvm(as.formula(txt.formula))
        lava::constrain(m,beta~alpha) <- function(x){-x}

        coef.alpha <- paste0(var.response,"~",var.breakpoint)
        coef.beta <- c(paste0(var.response,"~Us"),paste0(var.response,"~Us2"))    
        coef.gamma <- c(paste0(var.response,"~Vs"),paste0(var.response,"~Vs2"))
    }else if(pattern == "11"){
        formula.updated <- update(formula(object),".~.+Us+Vs")
        
        coef.alpha <- var.breakpoint
        coef.beta <- c("Us")
        coef.gamma <- c("Vs")
    }else if(pattern == "01"){
        formula.updated <- update(formula(object),".~Us+Vs")
        
        coef.alpha <- NULL
        coef.beta <- c("Us")
        coef.gamma <- c("Vs")
    }else if(pattern == "10"){
        txt.formula <- paste0(var.response,"~alpha*",var.breakpoint,"+beta*Us+Vs")
        m <- lava::lvm(as.formula(txt.formula))
        lava::constrain(m,beta~alpha) <- function(x){-x}

        coef.alpha <- paste0(var.response,"~",var.breakpoint)
        coef.beta <- paste0(var.response,"~Us")
        coef.gamma <- paste0(var.response,"~Vs")        
    }
    
    ## ** add variable
    iBreakpoint <- init.breakpoint
    if(trace){
        cat("Initialize breakpoints: ",paste(round(init.breakpoint, digits = digits), collapse = " "),"\n")
    }

    
    for(iIter in 1:n.iter){
        iBreakpointM1 <- iBreakpoint

        ## update design
        data[,c("Us") := (.SD[[1]]>iBreakpoint[1])*(.SD[[1]]-iBreakpoint[1]), .SDcols = var.breakpoint]
        data[,c("Vs") := -as.numeric(.SD[[1]]>iBreakpoint[1]), .SDcols = var.breakpoint]

        if(n.breakpoint>1){
            data[,c("Us2") := (.SD[[1]]>iBreakpoint[2])*(.SD[[1]]-iBreakpoint[2]), .SDcols = var.breakpoint]
            data[,c("Vs2") := -as.numeric(.SD[[1]]>iBreakpoint[2]), .SDcols = var.breakpoint]
        }
        
        ## *** estimate model coefficients
        if(pattern %in% c("111","11","01")){
            ebp <- lm(formula.updated, data = data)
        }else{
            ebp <- lava::estimate(m, data = data)
        }

        ## update breakpoint
        coef.tempo <- summary(ebp)$coef[,"Estimate"]
        iBreakpoint <- coef.tempo[coef.gamma]/coef.tempo[coef.beta] + iBreakpoint

        if(trace){
            cat("Iteration ",iIter,", breakpoints: ",paste(round(iBreakpoint, digits = digits), collapse = " "),"\n", sep = "")
        }

        ## exit
        iDiff <- abs(iBreakpoint-iBreakpointM1)
        if(any(iDiff<tol)){
            break
        }
        
    }

    ## ** standard error
    se <- rep(NA, n.breakpoint)
    if(pattern %in% c("111","11","01")){        
        vcov.tempo <- vcov(ebp)
    }else if(pattern %in% c("101","10")){
        iid.tempo <- iid(ebp)

        vec.sigma2 <- diag(vcov(ebp))
        vec.sigma2.robust <- colSums(iid.tempo^2)
        iid.tempo.scaled <- sweep(iid.tempo, FUN = "*", MARGIN = 2, STATS = sqrt(vec.sigma2/vec.sigma2.robust))
        colnames(iid.tempo.scaled)
        col.tempo <- -iid.tempo.scaled[,coef.alpha[1],drop =FALSE]
        colnames(col.tempo) <- coef.beta[1]

        vcov.tempo <- crossprod(cbind(iid.tempo.scaled, col.tempo))

        ## SE via the influence function
        ## term1 <- iid.tempo.scaled[,coef.gamma[1]]/coef(ebp)[coef.beta[1]] 
        ## term2 <- - coef(ebp)[coef.gamma[1]]* iid.tempo.scaled[,coef.beta[1]]/coef(ebp)[coef.beta[1]]^2
        ## sqrt(sum((term1 + term2)^2))
    }
    
    for(iBP in 1:n.breakpoint){
        term1 <- vcov.tempo[coef.gamma[iBP],coef.gamma[iBP]] / coef.tempo[coef.beta[iBP]]^2
        term2 <- vcov.tempo[coef.beta[iBP],coef.beta[iBP]] * (coef.tempo[coef.gamma[iBP]]/coef.tempo[coef.beta[iBP]]^2)^2
        term3 <- -2*vcov.tempo[coef.gamma[iBP],coef.beta[iBP]] * coef.tempo[coef.gamma[iBP]]/coef.tempo[coef.beta[iBP]]^3
        se[iBP] <- sqrt(term1+term2+term3)
    }
    
    ## ** fitted
    data.fit <- data.table(seq(min(data[[var.breakpoint]]), max(data[[var.breakpoint]]), length.out = n.points))
    names(data.fit) <- var.breakpoint

    data.fit[,c("Us") := (.SD[[1]]>iBreakpoint[1])*(.SD[[1]]-iBreakpoint[1]), .SDcols = var.breakpoint]
    data.fit[,c("Vs") := -as.numeric(.SD[[1]]>iBreakpoint[1]), .SDcols = var.breakpoint]
    if(n.breakpoint == 2){            
        data.fit[,c("Us2") := (.SD[[1]]>iBreakpoint[2])*(.SD[[1]]-iBreakpoint[2]), .SDcols = var.breakpoint]
        data.fit[,c("Vs2") := -as.numeric(.SD[[1]]>iBreakpoint[2]), .SDcols = var.breakpoint]
    }

    if(pattern %in% c("111","11","01")){
        data.fit[,fit := predict(ebp, newdata = data.fit)]
    }else{
        data.fit[,fit := predict(ebp, data = data.fit)[,1]] ## not newdata because it is lava::predict.lvmfit
    }
    
    ## ** display
    fit <- NULL ## [:CRANtest:] ggplot2
    gg <- ggplot2::ggplot(mapping = aes_string(var.breakpoint))
    gg <- gg + ggplot2::geom_point(data = cbind(data, observation = "observation"), aes_string(y = var.response, color = "observation"))
    gg <- gg + ggplot2::geom_line(data = data.fit, aes(y = fit, color = "fit"))
    gg <- gg + ggplot2::scale_colour_manual(name = "",
                                   values = c("red","black"))
    gg <- gg + ggplot2::ggtitle(label = paste0("pattern: ",pattern," | residual difference: ",round(iDiff, digits = digits)))

    ## ** export
    breakpoint <- list(value = iBreakpoint,
                       se = se,
                       init = init.breakpoint,
                       fit = data.fit[],
                       plot = gg,
                       opt <- c("n.iteration" = iIter,
                                "n.iterationMax" = n.iter,
                                "diff" = iDiff,
                                "tolerance" = tol)
                       )

    return(breakpoint)

}

######################################################################
### breakpoint.R ends here
