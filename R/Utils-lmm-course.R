### Utils-lmm-course.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  6 2019 (15:49) 
## Version: 
## Last-Updated: dec  6 2019 (18:18) 
##           By: Brice Ozenne
##     Update #: 114
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * procSummary
##' @title Compute summary statistics
##' @description Compute summary statistics (similar to the SAS macro procmean).
##' This is essentially an interface to the \code{stats::aggregate} function.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender + Gene will compute for each gender and gene the summary statistics for Y1 and for Y2.
##' Passed to the \code{stats::aggregate} function.
##' @param data [data.frame] dataset (in the wide format) containing the observations.
##' @param na.action [function] a function which indicates what should happen when the data contain 'NA' values.
##' Passed to the \code{stats::aggregate} function.
##' @param na.rm [logical] Should the summary statistics be computed by omitting the missing values.
##' @param which [character vector] name of the summary statistics to kept in the output.
##' Can be any of, or a combination of: \code{"observed"} (number of observations with a measurement),
##' \code{"missing"} (number of observations with a missing value), \code{"mean"}, \code{"sd"}, \code{"min"}, \code{"median"}, \code{"max"}.
##'
##' @examples
##' ## simulate data in the wide format
##' library(lava)
##' m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
##' categorical(m, labels = c("male","female")) <- ~gender
##' transform(m, id~gender) <- function(x){1:NROW(x)}
##' distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)
##'
##' set.seed(10)
##' d <- lava::sim(m, 1e2)
##'
##' ## add a missing value
##' d2 <- d
##' d2[1,"Y2"] <- NA
##'
##' ## run procSummary
##' procSummary(Y1+Y2 ~ 1, data = d)
##' procSummary(Y1+Y2 ~ gender, data = d)
##' 
##' procSummary(Y1+Y2 ~ gender, data = d2)
##' procSummary(Y1+Y2 ~ gender, data = d2, na.rm = TRUE)
##' @export
procSummary <- function(formula, data, na.action = na.pass, na.rm = FALSE,
                     which = c("observed","missing","mean","sd","min","median","max")){

    require(formula.tools)
    
    ## ** check and normalize user imput
    valid.which <- c("observed","missing","mean","sd","min","median","max")
    
    which <- match.arg(which, choices = valid.which, several.ok = TRUE)
    name.all <- all.vars(formula)
    if(any(name.all %in% c(valid.which, "outcome"))){
        invalid <- name.all[name.all %in% c(valid.which, "outcome")]
        stop("Name(s) \"",paste(invalid, collapse = "\" \""),"\" are used internally. \n",
             "Consider renaming the variables in the dataset and updating the formula. \n",
             sep = "")
    }
    if(any(name.all %in% names(data) == FALSE)){
        invalid <- name.all[name.all %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }
    
    name.Y <- formula.tools::lhs.vars(formula)
    n.Y <- length(name.Y)
    if(n.Y==0){
        stop("Wrong specification of argument \'formula\'. \n",
             "there need to be at least one variable in the left hand side of the formula. \n")
    }
    
    name.X <- formula.tools::rhs.vars(formula)
    n.X <- length(name.Y)

    
    ## ** compute summary statistics
    out <- NULL
    for(iY in 1:n.Y){ ## iY <- 1
        iFormula <- update(formula, paste0(name.Y[iY],"~."))
        iAggregate <- aggregate(iFormula, data=data, function(x){
            c("observed" = sum(!is.na(x)),
              "missing" = sum(is.na(x)),
              "mean" = mean(x, na.rm = na.rm),
              "sd" = sd(x, na.rm = na.rm),
              "min" = min(x, na.rm = na.rm),
              "median" = median(x, na.rm = na.rm),
              "max" = max(x, na.rm = na.rm))},
            na.action=na.action)
        iDF <- cbind(outcome = name.Y[iY], iAggregate[name.X], iAggregate[[name.Y[iY]]][,which,drop=FALSE])
        out <- rbind(out,iDF)
    }

    if(na.rm){
        attr(out,"correlation") <- cor(data[,name.Y], use = "pairwise")
    }else{
        attr(out,"correlation") <- cor(data[,name.Y])
    }
    
    ## ** export
    return(out)
}


## * lmm
##' @title Linear mixed model
##' @description Fit a linear mixed model using either a compound symmetry structure or an unstructured covariance matrix.
##' This is essentially an interface to the \code{nlme::gls} function.
##'
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param covariance [formula] Specify the model for the covariance.
##' No left hand side. On the right hand side,
##' either only the grouping variable (when specifying a compound symmetry structure), e.g. ~1|id,
##' or the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##'  
##' @examples
##' ## simulate data in the wide format
##' library(lava)
##' m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
##' categorical(m, labels = c("male","female")) <- ~gender
##' transform(m, id~gender) <- function(x){1:NROW(x)}
##' distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)
##'
##' set.seed(10)
##' dW <- lava::sim(m, 1e2)
##'
##' ## move to the long format
##' name.varying <- paste0("Y",1:4)
##' dL <- reshape(dW, direction  = "long",
##'               idvar = c("id","age","gender"),
##'               varying = name.varying,
##'               v.names = "Y",
##'               timevar = "visit")
##' rownames(dL) <- NULL
##' dL$visit <- factor(dL$visit,
##'                    levels = 1:length(name.varying),
##'                    labels = name.varying)
##' head(dL)
##' 
##' ## fit mixed model
##' e0.lmm <- lmm(Y ~ visit + age + gender, covariance = ~1|id, data = dL)
##' summary(e0.lmm)
##' nlme:::summary.gls(e0.lmm)
##'
##' e.lmm <- lmm(Y ~ visit + age + gender, covariance = ~visit|id, data = dL)
##' cat(attr(e.lmm,"code")) ## code used to fit the model
##' head(attr(e.lmm,"data")) ## data used to fit the model
##' summary(e.lmm)
lmm <- function(formula, covariance, data, df = FALSE, ...){

    ## ** check and normalize user imput
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' must be of class formula \n",
             "Something like: outcome ~ fixedEffect1 + fixedEffect2 \n")
    }
    name.mean <- all.vars(formula)
    if(any(name.mean %in% names(data) == FALSE)){
        invalid <- name.mean[name.mean %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    if(!inherits(covariance,"formula")){
        stop("Argument \'covariance\' must be of class formula. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    if(length(formula.tools::lhs.vars(covariance))!=0){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should not have any variable on the left hand side. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    if(!grepl("|",deparse(covariance),fixed = TRUE)){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "No | symbol found so no grouping variable could be defined. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    if(length(grepl("|",deparse(covariance),fixed = TRUE))>1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "The symbol | should only appear once. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    res.split <- strsplit(deparse(covariance),"|", fixed = TRUE)[[1]]
    name.cluster <- trimws(res.split[2], which = "both")
    name.repetition <- all.vars(as.formula(res.split[1]))
    if(length(name.repetition)>1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should have at most one variable before the grouping symbol (|). \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }

    ## ** fit mixed model
    if(length(name.repetition)==0){
        form.cor <- as.formula(paste0("~1|",name.cluster))
        
        e.lmm <- eval(parse(text = paste0("nlme::gls(",deparse(formula),",
                     correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                     data = data,
                     ...)")))
    }else{
        name.repetition.index <- paste0(name.repetition,".index")
        if(name.repetition.index %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable ",name.repetition.index," is used internally but already exists in \'data\' \n")
        }
        data[[name.repetition.index]] <- as.numeric(as.factor(data[[name.repetition]]))

        form.cor <- as.formula(paste0("~",name.repetition.index,"|",name.cluster))
        form.var <- as.formula(paste0("~1|",name.repetition))
        e.lmm <- eval(parse(text = paste0("nlme::gls(",deparse(formula),",
                     correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                     weights = nlme::varIdent(form = ",deparse(form.var),"),
                     data = data,
                     ...)")))
    }

    ## ** small sample correction
    if(df){
        stop("not implemented yet!")
        require(lavaSearch2)
        sCorrect(e.lmm) <- TRUE
    }

    ## ** export
    attr(e.lmm, "code") <- paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(e.lmm$call), collapse = ""))),"\n")
    attr(e.lmm, "data") <- data
    class(e.lmm) <- append("lmm",class(e.lmm))
    return(e.lmm)
}

## * summary.lmm
##' @title Summary Output for a Linear Mixed Model 
##' @description Summary output for a linear mixed model fitted with \code{lmm}.
##' This is a modified version of the \code{nlme::summary.gls} function.
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param digit [integer,>0] number of digit used to display numeric values.
##' @param conf.level [numeric,0-1] confidence level for the confidence intervals. 
summary.lmm <- function(object, digit = 3, conf.level = 0.95){

    ## ** welcome message
    if(!is.null(object$modelStruct$varStruct)){
        cat("  Linear mixed effect model with an unstructured covariance matrix \n")
    }else{
        cat("  Linear mixed effect model with a compound symmetry covariance matrix \n")
    }
    if(object$dim$REML){
        cat("  - fitted using Restricted Maximum Likelihood (REML) \n")
    }else{
        cat("  - fitted using Maximum Likelihood (REML) \n")
    }
    cat("  - likelihood :", as.double(object$logLik), " (df = ",object$dim$p,")\n",sep="")    
    cat(" \n")

    cat("Dataset:", deparse(object$call$data), "\n")
    cat(" - ", attr(object$modelStruct$corStruct,"Dim")$M, " clusters \n" , sep = "")
    cat(" - ", attr(object$modelStruct$corStruct,"Dim")$N, " observations \n",  sep = "")
    cat(" - ", attr(object$modelStruct$corStruct,"Dim")$maxLen, " maximum number of observations per cluster \n", sep = "")
    if(length(object$contrasts)>0){
        cat(" - levels of the categorical variables \n", sep = "")
        print(object$contrasts)
        
        if(attr(terms(eval(object$call$model)),"intercept") == 1){
            ref.level <- paste(unlist(lapply(names(object$contrasts), function(iC){
                paste0(iC,"=",rownames(object$contrasts[[iC]])[1])
            })), collapse = " ; ")
            cat(" - reference level: ",ref.level," \n", sep = "")
            cat(" \n")
        }
    }else{
        cat(" \n")

    }

    ## ** correlation structure
    cat("Correlation structure:",deparse(object$call$correlation),"\n")
    vec.corcoef <- coef(object$modelStruct$corStruct, unconstrained = FALSE)

    n.rep <- attr(object$modelStruct$corStruct,"Dim")$maxLen
    if(!is.null(object$modelStruct$varStruct)){
        name.rep <- attr(object$modelStruct$varStruct,"groupNames")
    }else{
        name.rep <- NULL
    }
    M.corcoef <- matrix(NA, nrow = n.rep, ncol = n.rep,
                        dimnames = list(name.rep,name.rep))
    M.corcoef[lower.tri(M.corcoef)] <- vec.corcoef
    M.corcoef[upper.tri(M.corcoef)] <- t(M.corcoef)[upper.tri(M.corcoef)]
    diag(M.corcoef) <- 1
    print(M.corcoef, digit = digit)
    cat("\n")

    ## ** variance structure
    cat("Variance structure:",deparse(object$call$weights),"\n")
    M.varcoef <- matrix(NA, nrow = 2, ncol = n.rep,
                        dimnames = list(c("variance","relative variance"),name.rep))
    if(!is.null(object$modelStruct$varStruct)){
        vec.varcoef <- coef(object$modelStruct$varStruct, unconstrained = FALSE)^2
        M.varcoef[1,] <- sigma(object)^2*c(1,vec.varcoef)
        M.varcoef[2,] <- c(1,vec.varcoef)
    }else{
        M.varcoef[1,] <- sigma(object)^2
        M.varcoef[2,] <- 1
    }
    print(M.varcoef, digit = digit)
    cat("\n")

    ## ** mean structure
    cat("Mean structure:",deparse(object$call$model),"\n")
    object2 <- object
    class(object2) <- setdiff(class(object2),"lmm")
    tTable <- summary(object2, verbose=FALSE)$tTable
    starSymbol <- symnum(tTable[,4], corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
    tTable <- data.frame(tTable[,1], nlme::intervals(object, which = "coef", level = conf.level)$coef[,c("lower","upper")],tTable[,2:4],starSymbol)
    names(tTable) <- c("estimate","lower","upper","se","t-value","p-value","")
    tTable.print <- tTable[,-5]
    tTable.print[["p-value"]] <- format.pval(tTable[["p-value"]], digits = digit, eps = 10^(-digit))
    tTable.print[["estimate"]] <- as.character(round(tTable[["estimate"]], digits = digit))
    tTable.print[["lower"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["upper"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["se"]] <- as.character(round(tTable[["se"]], digits = digit))
    print(tTable.print)
    cat("\n")
    cat("The lower and upper columns corresponds to the ",100*conf.level,"% confidence interval of the estimated coefficient\n", sep = "")
    cat("Note: p-value(s) and confidence interval(s) are not adjusted for multiple comparisons\n")

    ## ** export
    return(invisible(list(correlation = M.corcoef,
                          variance = M.varcoef,
                          mean = tTable)))
}

######################################################################
### Utils-lmm-course.R ends here
