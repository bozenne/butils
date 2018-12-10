### massLMorLVM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 10 2018 (16:12) 
## Version: 
## Last-Updated: dec 10 2018 (17:13) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Multiple Linear Regressions or Single LVM
#' @description Multiple Linear Regressions or Single LVM
#'
#' @param formula [formula or list of formula] formula whose right hand side indicates the covariantes.
#' The left hand side is ignored.
#' @param data [data.frame] dataset
#' @param exposure [character] the name of one variable in the formula whose effect should be assessed.
#' @param outcomes [character vector] name of the outcomes.
#' @param ... arguments passed to \code{glht2}.
#' 
#' @examples
#' \dontrun{
#' mSim <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ 2*eta + age, eta ~ E)
#' categorical(mSim, labels = c("yes","no")) <- ~E
#' latent(mSim) <- ~eta
#'
#' d <- lava::sim(mSim, n = 1e2, latent = FALSE)
#'
#' out <- massLMorLVM(~age+E, data = d, exposure = "E",
#'                    outcomes = c("Y1","Y2","Y3","Y4","Y5"))
#'
#' summary(out$model$lm[[1]])$coefficient["Eno",,drop=FALSE]
#' summary(out$model$lvm)$coef["Y1~Eno",,drop=FALSE]
#' c(var.lm = sigma(out$model$lm[[1]])^2,
#'   var.lvm = as.double(coef(out$model$lvm)["Y1~~Y1"] + coef(out$model$lvm)["eta~~eta"]))
#' 
#' lapply(out$glht, summary)
#'
#' dL <- melt(cbind(id = 1:NROW(d),d), id.vars = c("id","E","age"))
#' dL$variable <- as.factor(dL$variable)
#' e.gls <- gls(value ~ variable*E + variable*age,
#'              correlation = corSymm(form =~ as.numeric(variable)|id),
#'              weights = varIdent(form =~ 1|variable),
#'              data = dL)
#' summary(e.gls)$tTable
#' }
#' 
massLMorLVM <- function(formula, data, exposure, outcomes = NULL,
                        ...){

    ## ** normalize parameters
    if(inherits(formula, "formula")){

        if(is.null(outcomes)){
            stop("Argument \'outcomes\' must be specified when argument \'formula\' is a formula")
        }        
        formula <- lapply(outcomes, function(iY){
            update(formula,paste0(iY,"~."))
        })
        names(formula) <- outcomes
        
    }else if(inherits(formula, "list")){

        test.formula <- unlist(lapply(formula, function(x){inherits(x,"formula")}))
        if(any(test.formula == FALSE)){
            stop("Argument \'formula\' must be a formula or a list of formula \n")
        }else if(is.null(outcomes)){
            outcomes <- unlist(lapply(formula, function(iY){all.vars(iY)[1]}))
        }else if(length(outcomes) != length(formula)){
            stop("Argument \'formula\' and argument \'outcomes\' must have same length \n")
        }
        
    }else{
        stop("Argument \'formula\' must be a formula or a list of formula \n")
    }

    all.vars <- unique(unlist(lapply(formula, all.vars)))

    if(any(all.vars %in% names(data) == FALSE)){
        txt <- all.vars[all.vars %in% names(data) == FALSE]
        stop("Argument \'formula\' or \'outcomes\' inconsistent with \'data\' \n",
             "unknown variables: \"",paste(txt, collapse ="\" \""),"\"\n")
    }
    all.covars <- unique(unlist(lapply(formula, function(x){all.vars(update(x,"0~."))})))
    if(exposure %in% all.covars == FALSE){
        stop("Argument \'exposure\' is not among the variable defined by argument \'formula\' \n")
    }
    test.eta <- any(all.vars == "eta")
    if(any(test.eta)){
        stop("Argument \'formula\' must not contain any variable called eta \n")
    }
    
    n.outcomes <- length(outcomes)
    formula.lava <- lapply(formula, function(iF){update(iF,".~.+eta")})
    
    ## ** mass LM
    ls.lm <- lapply(formula, lm, data = data)
    class(ls.lm) <- "mmm"
    C.lm <- lavaSearch2::createContrast(ls.lm, var = exposure, add.variance = TRUE)$contrast
    e.glht_lm <- lavaSearch2::glht2(ls.lm, linfct = C.lm, ...)
    
    ## ** LVM
    m <- lvm(formula.lava)
    latent(m) <- ~eta
    e.lvm <- estimate(m, data = data)

    C.lvm <- lavaSearch2::createContrast(e.lvm, var = exposure)$contrast
    e.glht_lvm <- lavaSearch2::glht2(e.lvm, linfct = C.lvm, rhs = rep(0,NROW(C.lvm)), ...)


    ## ** output
    return(list(model = list(lm = ls.lm,
                             lvm = e.lvm),
                glht = list(lm = e.glht_lm,
                            lvm = e.glht_lvm)))
}

######################################################################
### massLMorLVM.R ends here
