### EmpLikRRCI.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 29 2018 (16:08) 
## Version: 
## Last-Updated: Apr  4 2023 (21:02) 
##           By: Paul Blanche
##     Update #: 81
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


EmpLikRRCI <- function (tstar,  data, alpha=0.05,
                        mytol=10^(-5),
                        maxupper=10,
                        minlower=0.1,
                        minupper=NULL,
                        maxlower=NULL,                        
                        whichAlgo=1
                        ){   
    # without constraint
    AJfit0 <- AJwithSE(tstar=tstar,data=data[data$group==0,])
    CIF1star0 <- AJfit0$res.tstar$risk.t
    AJfit1 <- AJwithSE(tstar=tstar,data=data[data$group==1,])
    CIF1star1 <- AJfit1$res.tstar$risk.t
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
    CIFupper <- NA
    try(CIFupper <- stats::uniroot(f=function(x){
        LikeRatioRRRdiff(data=data, RR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minupper,
        upper=maxupper, 
        tol=mytol)$root, silent = TRUE)
    # compute lower
    CIFlower <- NA
    try(CIFlower <- stats::uniroot(f=function(x){
        LikeRatioRRRdiff(data=data, RR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minlower,
        upper=maxlower,
        tol=mytol)$root, silent = TRUE)
    #--
    CI <- c(CIFlower,CIFupper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}

##----------------------------------------------------------------------
### EmpLikRRCI.R ends here
