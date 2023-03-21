### EmpLikSRCIKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec 20 2021 (14:22) 
## Version: 
## Last-Updated: Mar 21 2023 (14:28) 
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

EmpLikSRCIKM <- function (tstar,  data, alpha=0.05,
                          mytol=10^(-5),
                          maxupper=10,
                          minlower=0.1,
                          minupper,
                          maxlower,                        
                          whichAlgo=2
                          ){   
    # without constraint
    d0 <- data[data$group==0,]
    d1 <- data[data$group==1,]
    #--
    KMfit0 <- KMwithSE(tstar=tstar,data=d0)
    KMfit1 <- KMwithSE(tstar=tstar,data=d1)   
    Surv1star0 <- KMfit1$res.tstar$Surv
    Surv1star1 <- KMfit0$res.tstar$Surv
    #--
    unconstrSR <- Surv1star1/Surv1star0
    #---    
    if(missing(minupper)){
        minupper <- unconstrSR + mytol
    }
    if(missing(maxlower)){
        maxlower <- unconstrSR - mytol
    }
    #----
    # compute lower
    upper <- NA
    try(upper <- stats::uniroot(f=function(x){
        LikeRatioRRRdiffKM(data=data, SR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minupper,
        upper=maxupper, 
        tol=mytol)$root, silent = TRUE)
    # compute lower    
    lower <- NA
    try(lower <- stats::uniroot(f=function(x){
        LikeRatioRRRdiffKM(data=data, SR=x, tstar=tstar,whichAlgo=whichAlgo)$pval - alpha},
        lower=minlower,
        upper=maxlower,
        tol=mytol)$root, silent = TRUE)
    # return results
    CI <- c(lower,upper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}



#----------------------------------------------------------------------
### EmpLikSRCIKM.R ends here
