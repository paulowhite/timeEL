### EmpLikRdiffCIKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec 20 2021 (14:12) 
## Version: 
## Last-Updated: Mar 21 2023 (14:26) 
##           By: Paul Blanche
##     Update #: 20
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

EmpLikRdiffCIKM <- function (tstar,  data, alpha=0.05,
                             mytol=10^(-5),
                             maxupper=0.99,
                             minlower=-0.99,
                             minupper,
                             maxlower,                        
                             whichAlgo=1
                             ){
    # without constraint
    KMfit0 <- KMwithSE(tstar=tstar,data=data[data$group==0,])
    CIF1star0 <- 1-KMfit0$res.tstar$Surv
    KMfit1 <- KMwithSE(tstar=tstar,data=data[data$group==1,])
    CIF1star1 <- 1-KMfit1$res.tstar$Surv
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

        LikeRatioRRRdiffKM(data=data, Rdiff=x, tstar=tstar,
                           whichAlgo=whichAlgo)$pval - alpha},
        lower=minupper,
        upper=maxupper, 
        tol=mytol
        )$root, silent = TRUE)   
    # compute lower
    CIFlower <- NA
    try(CIFlower <- stats::uniroot(f=function(x){
        LikeRatioRRRdiffKM(data=data, Rdiff=x, tstar=tstar, 
                           whichAlgo=whichAlgo)$pval - alpha},
        lower=minlower,
        upper=maxlower,
        tol=mytol)$root, silent = TRUE)
    # return results
    CI <- c(CIFlower,CIFupper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}


#----------------------------------------------------------------------
### EmpLikRdiffCIKM.R ends here
