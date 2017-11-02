### onload.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 16 2017 (23:42) 
## Version: 
## last-updated: nov  2 2017 (10:06) 
##           By: Brice Ozenne
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## Functions not exported by boot
boot_isMatrix <- get("isMatrix", envir = asNamespace("boot"), inherits = FALSE)
boot_index.array <- get("index.array", envir = asNamespace("boot"), inherits = FALSE)
boot_boot.return <- get("boot.return", envir = asNamespace("boot"), inherits = FALSE)

riskRegression_coxVariableName <- get("coxVariableName", envir = asNamespace("riskRegression"), inherits = FALSE)
riskRegression_coxDesign <- get("coxDesign", envir = asNamespace("riskRegression"), inherits = FALSE)


#----------------------------------------------------------------------
### onload.R ends here
