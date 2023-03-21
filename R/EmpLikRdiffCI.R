### EmpLikRdiffCI.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 29 2018 (16:08) 
## Version: 
## Last-Updated: Mar 21 2023 (14:28) 
##           By: Paul Blanche
##     Update #: 86
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


EmpLikRdiffCI <- function (tstar,  data, alpha=0.05,
                           mytol=10^(-5),
                           maxupper=0.99,
                           minlower=-0.99,
                           minupper,
                           maxlower,                        
                           whichAlgo=1
                           ){
    # without constraint
    res.AJwithSE.0 <- AJwithSE(tstar=tstar,data=data[data$group==0,],CompAllse=FALSE)
    res.AJwithSE.1 <- AJwithSE(tstar=tstar,data=data[data$group==1,],CompAllse=FALSE) 
    CIF1star0 <- res.AJwithSE.0$risk.t
    CIF1star1 <- res.AJwithSE.1$risk.t
    unconstrRdiff <- CIF1star1-CIF1star0    
    #---
    if(missing(minupper)){
        minupper <- unconstrRdiff + mytol
    }
    if(missing(maxlower)){
        maxlower <- unconstrRdiff - mytol
    }
    #----
    # compute lower
    CIFupper <- NA
    try(CIFupper <- stats::uniroot(f=function(x){
        LikeRatioRRRdiff(data=data, Rdiff=x, tstar=tstar,
                         whichAlgo=whichAlgo)$pval - alpha},
        lower=minupper,
        upper=maxupper, 
        tol=mytol
        )$root, silent = TRUE)   
    # compute lower
    CIFlower <- NA
    try(CIFlower <- stats::uniroot(f=function(x){
        LikeRatioRRRdiff(data=data, Rdiff=x, tstar=tstar, 
                         whichAlgo=whichAlgo)$pval - alpha},
        lower=minlower,
        upper=maxlower,
        tol=mytol)$root, silent = TRUE)
    #--
    CI <- c(CIFlower,CIFupper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}

##----------------------------------------------------------------------
### EmpLikRdiffCI.R ends here
