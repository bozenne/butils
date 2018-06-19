### onload.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 16 2017 (23:42) 
## Version: 
## last-updated: jun 19 2018 (10:22) 
##           By: Brice Ozenne
##     Update #: 8
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

estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE) ## needed for bootReg.lvm since estimate.lvm is saved in the call

riskRegression_coxVariableName <- get("coxVariableName", envir = asNamespace("riskRegression"), inherits = FALSE)
riskRegression_coxModelFrame <- get("coxModelFrame", envir = asNamespace("riskRegression"), inherits = FALSE)


#----------------------------------------------------------------------
### onload.R ends here
