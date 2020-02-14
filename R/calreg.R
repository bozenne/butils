### calreg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 14 2020 (17:23) 
## Version: 
## Last-Updated: feb 14 2020 (18:25) 
##           By: Brice Ozenne
##     Update #: 18
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * calreg (documentation)
#' @title Regression On Unobserved Exposure Using Calibration.
#' @description Perform a linear regression on an unobserved exposure (X) using a proxy (Z) whose relationship with the exposure has been studied using an external dataset.
#'
#' @param formula formula for the linear model.
#' @param data dataset used to fit the linear model relating the outcome (Y) to the (unobserved) exposure (X).
#' @param calibration a \code{lm} object or \code{nls} object relating a proxy (Z) to the unobserved exposure (X).
#'
#' @details Consider a first sample \eqn{(X_i,Z_i)_{i \in \{1,\ldots,m\}}} that is used to estimate \eqn{\alpha} in:
#' \deqn{X = f(alpha,Z) + \varepsilon_{\alpha}}
#' This is the model to give to the argument \code{calibration}.
#'
#' The aim is to use a second sample \eqn{(Y_j,Z_j)_{j \in \{1,\ldots,n\}}} to estimate \eqn{\beta_1} in:
#' \deqn{Y = \beta_0 + \beta_1 X + \varepsilon_{\beta}}
#' The formula of this model should be given to the argument \code{formula} and the dataset to the argument \code{data}.
#' This function uses the predictions of the first model to estimate \(\X\) for the second sample.
#' The uncertainty is decomposed into two part:
#' \itemize{
#' \item one related to the finite number of observations in the second sample.
#' \item one related to the estimation of the parameters in the calibration model, to account for the fact that \(X\) is estimated and not observed.
#' This part of the variance is estimated using a delta method.
#' }
#' 
#' @return A data.frame containing the estimates, standard errors, confidence intervals and p-values for each regression coefficient.
#' The output has an attribute \code{"regression"} containing the fitted linear model (ignoring the uncertainty related to the calibration)
#' and an attribute \code{"var.add"} representing additional variance-covariance matrix due to the calibration.
#'
#' @examples
#' library(lava)
#' n <- 1e2
#' 
#' ## linear case
#' mSim.lin <- lvm(fMRI ~ occ, occ ~ blood)
#' distribution(mSim.lin, ~blood) <- uniform.lvm(-0.9,2)
#'
#' set.seed(10)
#' d1.lin <- sim(mSim.lin, n = n)[,c("occ","blood"),drop=FALSE]
#' d2.lin <- sim(mSim.lin, n = n)[,c("fMRI","blood"),drop=FALSE]
#'
#' e1.lin <- lm(occ ~ blood, data = d1.lin)
#' res.lin <- calreg(fMRI ~ occ, data = d2.lin, calibration = e1.lin)
#' res.lin
#' summary(attr(res.lin, "regression"))$coef
#'
#' ## non-linear case
#' mSim.nlin <- lvm(fMRI ~ occ, occ[mu:0.1] ~ 0*blood)
#' distribution(mSim.nlin, ~blood) <- uniform.lvm(-0.9,3)
#' constrain(mSim.nlin, mu~blood) <- function(x){2*x/(1+x)}
#'
#' set.seed(10)
#' d1.nlin <- sim(mSim.nlin, n = n)[,c("occ","blood"),drop=FALSE]
#' d2.nlin <- sim(mSim.nlin, n = n)[,c("fMRI","blood"),drop=FALSE]
#'
#' ## gg <- ggplot(d1.nlin, aes(x = blood)) + geom_point(aes(y = occ))
#'
#' e1.nlin <- nls(occ ~ (occmax * blood)/(EC + blood), data = d1.nlin,
#'          start = list(occmax = 1, EC = 1))
#' d1.nlin$fit <- fitted(e1.nlin)
#' 
#' ## gg + geom_line(data = d1.nlin, aes(y = fit), color = "red")

## * calreg (code)
##' @rdname calreg
##' @export
calreg <- function(formula, data, calibration){

    ## ** check arguments
    if(!inherits(calibration,"lm") && !inherits(calibration,"nls")){
        stop("Only compatible with lm and nls objects \n")
    }

    ## ** extract information
    name.X <- all.vars(formula(calibration))[1]
    X <- model.matrix(calibration)
    alpha <- coef(calibration)
    vcov_alpha <- vcov(calibration)
    
    ## ** fit new regression
    refit <- function(p){
        data[[name.X]] <- X %*% p
        return(lm(formula, data = data))
    }
    regression <- refit(alpha)
    Sregression <- summary(regression)
    
    ## ** delta method
    dcoef <- numDeriv::jacobian(func = function(x){coef(refit(x))},
                                x = alpha)
    var.add <- dcoef %*% vcov_alpha %*% t(dcoef)
    ## dcoef[1,,drop=FALSE] %*% vcov_alpha %*% t(dcoef[1,,drop=FALSE]) - var.add[1,1]

    newcov <- Sregression$cov.unscaled * sigma(regression)^2 + var.add
    table2 <- data.frame(matrix(NA, nrow = nrow(Sregression$coef), ncol = 7,
                                dimnames = list(rownames(Sregression$coef), c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"))))
    table2$estimate <- as.double(Sregression$coef[,"Estimate"])
    table2$std.error <- as.double(sqrt(diag(newcov)))
    table2$df <- regression$df.residual
    table2$ci.lower <- table2$estimate + qt(0.025, df = table2$df) * table2$std.error
    table2$ci.upper <- table2$estimate + qt(0.975, df = table2$df) * table2$std.error
    table2$statistic <- table2$estimate/table2$std.error
    table2$p.value <- 2*(1-pt(abs(table2$statistic), df = table2$df))

    attr(table2, "regression") <- regression
    attr(table2, "var.add") <- var.add
    
    return(table2)
    
}

######################################################################
### calreg.R ends here
