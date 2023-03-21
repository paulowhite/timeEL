### EmpLikRRCIKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (13:02) 
## Version: 
## Last-Updated: Mar 21 2023 (14:27) 
##           By: Paul Blanche
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


EmpLikRRCIKM <- function (tstar,  data, alpha=0.05,
                          mytol=10^(-5),
                          maxupper=10,
                          minlower=0.1,
                          minupper,
                          maxlower,                        
                          whichAlgo=2
                          ){
    # without constraint
    KMfit0 <- KMwithSE(tstar=tstar,data=data[data$group==0,])
    CIF1star0 <- 1-KMfit0$res.tstar$Surv
    KMfit1 <- KMwithSE(tstar=tstar,data=data[data$group==1,])
    CIF1star1 <- 1-KMfit1$res.tstar$Surv
    unconstrRR <- CIF1star1/CIF1star0
    #---
    if(missing(minupper)){
        minupper <- unconstrRR + mytol
    }
    if(missing(maxlower)){
        maxlower <- unconstrRR - mytol
    }
    #----
    # compute lower
    upper <- NA
    try(upper <- stats::uniroot(f=function(x){
        LikeRatioRRRdiffKM(data=data, RR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minupper,
        upper=maxupper, 
        tol=mytol)$root, silent = TRUE)
    # compute lower    
    lower <- NA
    try(lower <- stats::uniroot(f=function(x){
        LikeRatioRRRdiffKM(data=data, RR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minlower,
        upper=maxlower,
        tol=mytol)$root, silent = TRUE)
    # return results
    CI <- c(lower,upper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}

#----------------------------------------------------------------------
### EmpLikRRCIKM.R ends here
