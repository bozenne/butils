### breakpoint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 26 2018 (09:13) 
## Version: 
## Last-Updated: jun 26 2018 (15:35) 
##           By: Brice Ozenne
##     Update #: 214
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
#' @param breakpoint.var [character] the variable regarding which the breakpoints should be found.
#' @param breakpoint.init [list] a list containing initialization values for the breakpoints for each pattern.
#' Initialisation values can be a vector or a matrix with the same number of rows as the number of breakpoints.
#' @param n.iter [integer, >0] the maximum number of iterations used to estimates the breakpoints.
#' @param n.init [integer, >0] the number of quantiles used to generate initialisation points.
#' Only active when breakpoint.init is \code{NULL}.
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
#' BIC(resBP)
#' gg + geom_line(data = resBP$BP101$fit, aes(y = fit))
#'
#' #### example from the package segmented
#' if(require(segmented)){
#' GS <- segmented(e.lm, psi = c(1,2))
#'
#' 
#' cbind(value = resBP$BP111$breakpoint,
#'       se = resBP$BP111$breakpoint.se)
#' GS$psi
#' }
#' 
#' if(require(gridExtra)){
#'   autoplot(resBP)
#' }
#'  

## * breakpoint (code)
#' @rdname breakpoint
#' @export
breakpoint <- function(object, pattern = c("111","101","10","11","01"), breakpoint.var = NULL, breakpoint.init = NULL,
                       n.iter = 10, tol = 1e-3, n.init = 5, n.points = 5,
                       trace = 2, digits = -log10(tol)){


    ## ** find data
    dt <- model.frame(object)
    if(!is.data.table(dt)){
        data.table::setDT(dt)
    }else{
        dt <- data.table::copy(dt)
    }

    ## ** check arguments
    ## *** object
    if("lm" %in% class(object) == FALSE){
        stop("Only work with lm objects \n")
    }
    
    ## *** data
    reserved.names <- c("fit","Us","Vs","Us2","Vs2","gamma","beta")
    if(any(names(dt) %in% reserved.names)){
        txt <- names(dt[names(dt) %in% reserved.names])
        stop("data contains reserved names: \"",paste0(txt, collapse = "\" \""),"\"\n")
    }


    ## *** response.var
    response.var <- all.vars(update(formula(object),".~0"))

    ## *** breakpoint.var
    if(is.null(breakpoint.var)){
        breakpoint.var <- all.vars(update(formula(object),"0~."))
    }

    if(length(breakpoint.var)!=1){
        stop("Argument \"breakpoint.var\" should correspond to only one variable \n")
    }

    if(breakpoint.var %in% names(dt) == FALSE){
        stop("Variable \"",breakpoint.var,"\" is not in data \n")
    }

    ## *** pattern
    if(any(pattern %in% c("111","101","10","11","01") == FALSE) ){
        stop("n.breakpoint must be one of \"111\",\"101\",\"10\",\"11\",\"01\" \n")
    }

    ## *** breakpoint.init
    if(is.null(breakpoint.init)){
        breakpoint.init <- list()
    }else{
        if(!is.list(breakpoint.init)){
            stop("Argument \'breakpoint.init\' must be a list \n")
        }
        if(any(pattern %in% names(breakpoint.init) == FALSE)){
            stop("Argument \'breakpoint.init\' must be a list containing elements \"",paste(pattern, collapse = "\" \""),"\" \n")
        }
    }

    ## ** compute breakpoints
    out <- vector(mode = "list", length = length(pattern))
    names(out) <- paste0("BP",pattern)
    for(iPattern in pattern){
        if(trace>0){ cat("* Search breaking points (pattern ",iPattern,")", sep = "") }
        if(trace>2){ cat("\n") }
        out[[paste0("BP",iPattern)]] <- .BPfit(object = object,
                                               data = dt,
                                               response.var = response.var,
                                               breakpoint.var = breakpoint.var,
                                               breakpoint.init = breakpoint.init[[iPattern]],
                                               pattern = iPattern,
                                               n.iter = n.iter,
                                               tol = tol,
                                               n.init = n.init,
                                               n.points = n.points,
                                               trace = (trace-1),
                                               digits = digits)
        if(trace <= 2){
            cat(" - cv=",out[[paste0("BP",iPattern)]]$cv," \n",sep="")
        }
        ## object$BP111$plot
    }
    
    ## ** export
    out$breakpoint.var <- breakpoint.var
    out$response.var <- response.var
    out$data <- dt
    class(out) <- "breakpoint"
    return(out)
}

## * BIC
#' @method BIC breakpoint
#' @export
BIC.breakpoint <- function(object,...){

    index.BP <- grep("^BP", names(object))
    
    return(unlist(lapply(object[index.BP],function(iP){
        if(!is.null(iP$model)){
            return(BIC(iP$model))
        }else{
            return(NA)
        }
    })))

}

## * autoplot
#' @title Display Regression Line and Observations
#' @description Display regression line and observations.
#'
#' @param object output of \code{breakpoint}
#' @param pattern [vector of character] the number and type of breakpoints to be display.
#' @param combine.plot [logical] should the plots for the different patterns be combined into one.
#' @param nrow [integer, >0] number of rows used when combining the plots.
#' @param ncol [integer, >0] number of columns used when combining the plots.
#' @param title [character] the title of the combined plot. 
#' @param plot [logical] should the plot be displayed in a window?
#' @param text.size [numeric, >0] the size of the text in the plot.
#' @param add.cv.title [logical] should the convergence status of the estimation algorithm
#' be displayed in the title of the plot.
#' @param add.bic.title [logical] should the bic of the model
#' be displayed in the title of the plot.
#' @param ... not used. For compatibility with the generic function.
#' 
#' 
#' @method autoplot breakpoint
#' @export 
autoplot.breakpoint <- function(object, pattern = NULL, plot = TRUE,
                                combine.plot = TRUE, nrow = NULL, ncol = NULL, title = NULL, text.size = 10,
                                add.cv.title = TRUE, add.bic.title = FALSE, ...){

    fit <- NULL ## [:CRANtest:] ggplot2

    ## ** normalize argument
    if(is.null(pattern)){
        pattern <- grep("^BP", names(object), value = TRUE)
    }else{
        pattern <- match.arg(pattern, choices = names(object), several.ok = TRUE)
    }
    n.pattern <- length(pattern)

    breakpoint.var <- object$breakpoint.var
    response.var <- object$response.var
    data <- object$data

    ##
    tryPkg <- requireNamespace("gridExtra")
    if("try-error" %in% class(tryPkg)){
        stop(tryPkg)
    }


    ## ** make individual plots
    ls.plot <- vector(mode = "list", length = n.pattern)
    names(ls.plot) <- pattern
    object.BIC <- BIC(object)
    
    for(iPattern in pattern){ ## iPattern <- pattern[2]
        title.txt <- paste0("pattern: ",iPattern)

        if(add.cv.title){title.txt <- paste0(title.txt," | convergence: ",object[[iPattern]]$cv)}
        if(add.bic.title){title.txt <- paste0(title.txt," | BIC: ",round(object.BIC[iPattern],3))}

        ls.plot[[iPattern]] <- ggplot2::ggplot(mapping = aes_string(breakpoint.var))
        ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::geom_point(data = cbind(data, observation = "observation"), aes_string(y = response.var, color = "observation"))
        if(!is.null(object[[iPattern]]$fit)){
            ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::geom_line(data = object[[iPattern]]$fit, aes(y = fit, color = "fit"))
            ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::scale_colour_manual(name = "",
                                                                                      values = c("red","black"))
        }else{
               ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::scale_colour_manual(name = "", values = c("black"))
        }
        
        ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::ggtitle(label = title.txt)
        ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::theme(text = element_text(size = text.size))

        if(pattern[which.min(object.BIC)]==iPattern){
            ls.plot[[iPattern]] <- ls.plot[[iPattern]] + ggplot2::theme(plot.title = element_text(colour = "darkblue"))
        }
    }

    if(combine.plot && n.pattern>1){
        vec.txt <- paste0("ls.plot[[\"",pattern,"\"]] + theme(legend.position=\"none\")")
        txt <- paste0(vec.txt, collapse = ", \n", sep = "")
        txt.all <- paste0("gridExtra::arrangeGrob(",txt,
                          if(!is.null(nrow)){paste0(", nrow = ",nrow)},
                          if(!is.null(ncol)){paste0(", ncol = ",ncol)},
                          if(!is.null(title)){paste0(",top=\"",title,"\"")},
                          ")")
        out <- eval(parse(text = txt.all))
        if(plot){
            gridExtra::grid.arrange(out)
        }
    }else{
        out <- ls.plot
        if(plot){
            lapply(out, print)
        }
    }

    
        return(invisible(out))
}

## * .calcBP
.BPfit <- function(object, data, response.var, breakpoint.var, breakpoint.init, pattern,
                   n.iter, tol, n.init, n.points,
                   trace, digits){

    n.breakpoint <- nchar(pattern)-1
    
    ## ** lvm model
    ## free flat free
    if(pattern == "111"){
        formula.updated <- update(formula(object),".~.+Us+Vs+Us2+Vs2")
        
        coef.beta <- c("Us","Us2")
        coef.gamma <- c("Vs","Vs2")
    }else if(pattern == "101"){
        formula.updated <- update(update(formula(object),paste0(".~.-",breakpoint.var)),
                                  paste0(".~.+I(",breakpoint.var,"-Us)+Vs+Us2+Vs2"))

        ## txt.formula <- paste0(response.var,"~alpha*",breakpoint.var,"+beta*Us+Vs+Us2+Vs2")
        ## m <- lava::lvm(as.formula(txt.formula))
        ## lava::constrain(m,beta~alpha) <- function(x){-x}

        ## coef.alpha <- paste0(response.var,"~",breakpoint.var)
        ## coef.beta <- c(paste0(response.var,"~Us"),paste0(response.var,"~Us2"))    
        ## coef.gamma <- c(paste0(response.var,"~Vs"),paste0(response.var,"~Vs2"))
        coef.beta <- c("Us","Us2")
        coef.gamma <- c("Vs","Vs2")

    }else if(pattern == "11"){
        formula.updated <- update(formula(object),".~.+Us+Vs")
        
        coef.beta <- c("Us")
        coef.gamma <- c("Vs")
    }else if(pattern == "01"){
        formula.updated <- update(update(formula(object),paste0(".~.-",breakpoint.var)),
                                  paste0(".~.+Us+Vs"))
        
        coef.beta <- c("Us")
        coef.gamma <- c("Vs")
    }else if(pattern == "10"){
        formula.updated <- update(update(formula(object),paste0(".~.-",breakpoint.var)),
                                  paste0(".~.+I(",breakpoint.var,"-Us)+Vs"))
        coef.beta <- c("Us")
        coef.gamma <- c("Vs")
    }
    
    ## ** loop over initializations points
    if(is.null(breakpoint.init)){
                
        probs.breakpoint <- seq(0,1, length.out = n.init+2)[2:(n.init+1)]
        quantile.data <- quantile(data[[breakpoint.var]], probs = probs.breakpoint)
        breakpoint.init <- utils::combn(quantile.data, m = n.breakpoint)

    }else {
        if(is.vector(breakpoint.init)){
            breakpoint.init <- cbind(as.double(breakpoint.init))
        }
        if(NROW(breakpoint.init) != n.breakpoint){
            stop("Incorrect initialization for the breakpoints (too many or to few) \n",
                 "NROW(breakpoint.init): ",NROW(breakpoint.init),"\n",
                 "required:",n.breakpoint,"\n")
        }
    }
    n.init <- NCOL(breakpoint.init)
    
    ls.res <- vector(mode = "list", length = n.init)
    vec.score <- rep(-Inf, times = n.init) 
    vec.score2 <- rep(-Inf, times = n.init)

    breakpoint.min <- min(data[[breakpoint.var]], na.rm = TRUE)
    breakpoint.max <- max(data[[breakpoint.var]], na.rm = TRUE)

    if(trace>0){cat(": ")}
    for(iInit in 1:n.init){
        ls.res[[iInit]] <- .warperBP(formula = formula.updated,
                                     data = data,
                                     breakpoint.init = breakpoint.init[,iInit],
                                     breakpoint.var = breakpoint.var,
                                     breakpoint.min = breakpoint.min,
                                     breakpoint.max = breakpoint.max,
                                     n.breakpoint = n.breakpoint,
                                     coef.beta = coef.beta,
                                     coef.gamma = coef.gamma,
                                     pattern = pattern,
                                     n.iter = n.iter,
                                     tol = tol,
                                     trace = (trace-1),
                                     digits = digits)

        if(ls.res[[iInit]]$cv){
            vec.score[iInit] <- logLik(ls.res[[iInit]]$model)
            if(ls.res[[iInit]]$continuity){
                vec.score2[iInit] <- logLik(ls.res[[iInit]]$model)
            }
        }
        if(trace>0){cat("+")}
    }

    ## ** find the best fit
    if(any(!is.infinite(vec.score2))){ ## if any model respect continuity constrain take the best of them
        out <- ls.res[[which.max(vec.score2)[1]]]
    }else{ ## other look over those with discontinuity
        out <- ls.res[[which.max(vec.score)[1]]]
    }
    out$all <- list(init = breakpoint.init,
                    score = vec.score,
                    breakpoint = do.call(rbind,lapply(ls.res,"[[","breakpoint"))
                    )
    
    ## ** standard error
    if(out$cv){
        out$breakpoint.se <- rep(NA, n.breakpoint)
        
        if(pattern %in% c("111","11","01")){        
            vcov.tempo <- vcov(out$model)
        }else if(pattern %in% c("101","10")){
            tryPkg <- requireNamespace("lavaSearch2")
            if("try-error" %in% class(tryPkg)){
                stop(tryPkg)
            }
            iid.tempo <- lavaSearch2::iid2(out$model, robust = FALSE)

            m.tempo <- cbind(-iid.tempo[,paste0("I(",breakpoint.var," - Us)")])
            colnames(m.tempo) <- "Us"
            vcov.tempo <- crossprod(cbind(iid.tempo,m.tempo))
        
        }
    
        ## SE via the influence function
        ## term1 <- iid.tempo.scaled[,coef.gamma[1]]/coef(ebp)[coef.beta[1]] 
        ## term2 <- - coef(ebp)[coef.gamma[1]]* iid.tempo.scaled[,coef.beta[1]]/coef(ebp)[coef.beta[1]]^2
        ## sqrt(sum((term1 + term2)^2))
    
        for(iBP in 1:n.breakpoint){
            term1 <- vcov.tempo[coef.gamma[iBP],coef.gamma[iBP]] / out$coef[coef.beta[iBP]]^2
            term2 <- vcov.tempo[coef.beta[iBP],coef.beta[iBP]] * (out$coef[coef.gamma[iBP]]/out$coef[coef.beta[iBP]]^2)^2
            term3 <- -2*vcov.tempo[coef.gamma[iBP],coef.beta[iBP]] * out$coef[coef.gamma[iBP]]/out$coef[coef.beta[iBP]]^3
            out$breakpoint.se[iBP] <- sqrt(term1+term2+term3)
        }
    }
    
    ## ** fitted
    if(!is.null(out$model)){

        if(n.breakpoint==1){
            out$fit <- data.table(c(seq(min(data[[breakpoint.var]]), out$breakpoint[1], length.out = n.points),
                                    seq(out$breakpoint[1], max(data[[breakpoint.var]]), length.out = n.points)))
        }else if(n.breakpoint==2){
            out$fit <- data.table(c(seq(min(data[[breakpoint.var]]), out$breakpoint[1], length.out = n.points),
                                    seq(out$breakpoint[1], out$breakpoint[2], length.out = n.points),
                                    seq(out$breakpoint[2], max(data[[breakpoint.var]]), length.out = n.points)))
        }
        names(out$fit) <- breakpoint.var

        out$fit[,c("Us") := (.SD[[1]]>out$breakpoint[1])*(.SD[[1]]-out$breakpoint[1]), .SDcols = breakpoint.var]
        out$fit[,c("Vs") := -as.numeric(.SD[[1]]>out$breakpoint[1]), .SDcols = breakpoint.var]
        if(n.breakpoint == 2){            
            out$fit[,c("Us2") := (.SD[[1]]>out$breakpoint[2])*(.SD[[1]]-out$breakpoint[2]), .SDcols = breakpoint.var]
            out$fit[,c("Vs2") := -as.numeric(.SD[[1]]>out$breakpoint[2]), .SDcols = breakpoint.var]
        }

        out$fit[,c("fit") := predict(out$model, newdata = out$fit)]
    }

    ## ** export
    return(out)

}

## * warperBP
.warperBP <- function(formula, data,
                      breakpoint.init, breakpoint.var, breakpoint.min, breakpoint.max, n.breakpoint,
                      coef.beta, coef.gamma, pattern,
                      n.iter, tol,
                      trace, digits){


    if(trace>0){
        cat("    Initialize breakpoints: ",paste(round(breakpoint.init, digits = digits), collapse = " "),"\n")
    }
    iBreakpoint <- breakpoint.init
    cv <- FALSE
    
    for(iIter in 1:n.iter){
        iBreakpointM1 <- iBreakpoint

        ## ** update design
        data[,c("Us") := (.SD[[1]]>iBreakpoint[1])*(.SD[[1]]-iBreakpoint[1]), .SDcols = breakpoint.var]
        data[,c("Vs") := -as.numeric(.SD[[1]]>iBreakpoint[1]), .SDcols = breakpoint.var]

        if(n.breakpoint>1){
            data[,c("Us2") := (.SD[[1]]>iBreakpoint[2])*(.SD[[1]]-iBreakpoint[2]), .SDcols = breakpoint.var]
            data[,c("Vs2") := -as.numeric(.SD[[1]]>iBreakpoint[2]), .SDcols = breakpoint.var]
        }
        
        ## ** estimate model coefficients
        ebp <- lm(formula, data = data)
        
        ## ** update breakpoint
        coef.tempo <- summary(ebp)$coef[,"Estimate"]
        if(pattern=="101"){
            coef.tempo["Us"] <- -coef.tempo[paste0("I(",breakpoint.var," - Us)")]
        }else if(pattern == "10"){
            coef.tempo["Us"] <- -coef.tempo[paste0("I(",breakpoint.var," - Us)")]
        }
        iBreakpoint <- coef.tempo[coef.gamma]/coef.tempo[coef.beta] + iBreakpoint
        
        if(any(is.na(iBreakpoint)) || any(iBreakpoint<breakpoint.min) || any(iBreakpoint>breakpoint.max) || is.unsorted(iBreakpoint) ){
            ## case where one breakpoint is outside the domain
            ## or when the breakpoints are not in increasing order
            return(list(model = NULL,
                        coef = coef.tempo,
                        breakpoint = iBreakpoint,
                        breakpoint.init = breakpoint.init,
                        cv = FALSE,
                        continuity = FALSE))
        }
        
        ## ** display
        if(trace>0){
            cat("    > iteration ",iIter,
                ", breakpoints: ",paste(round(iBreakpoint, digits = digits), collapse = " "),
                ", Vs: ",paste(round(coef.tempo[coef.gamma], digits = digits), collapse = " "),
                " \n", sep = "")
        }

        ## ** cv
        iDiff <- abs(iBreakpoint-iBreakpointM1)
        if(all(iDiff<tol)){
            cv <- TRUE
            break
        }
        
    }

    ## ** enforce continuity
    test.continuity <- all(abs(coef.tempo[coef.gamma])<tol)
    if(test.continuity==FALSE){
        ebp <- lm(update(formula,".~.-Vs-Vs2"), data = data)

        coef.tempo <- summary(ebp)$coef[,"Estimate"]
        if(pattern=="101"){
            coef.tempo["Us"] <- -coef.tempo[paste0("I(",breakpoint.var," - Us)")]
        }else if(pattern == "10"){
            coef.tempo["Us"] <- -coef.tempo[paste0("I(",breakpoint.var," - Us)")]
        }

    }
    
    ## ** export
    test.continuity <- all(abs(coef.tempo[coef.gamma])<tol)
    
    return(list(model = ebp,
                coef = coef.tempo,
                breakpoint = iBreakpoint,
                breakpoint.init = breakpoint.init,
                diff = iDiff,
                n.iter = n.iter,
                cv = cv,
                continuity = test.continuity))
}

######################################################################
### breakpoint.R ends here
