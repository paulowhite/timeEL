### LRCIKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (10:27) 
## Version: 
## Last-Updated: Mar 21 2023 (14:27) 
##           By: Paul Blanche
##     Update #: 78
#----------------------------------------------------------------------
## 
### Commentary: 
##
##  Empirical likelihood (ratio) based confidence
##  interval for usual survival data (Kaplan-Meier, no competing risks)
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


LRCIKM <- function (tstar,
                    data,
                    alpha=0.05,
                    mytol=10^(-4),
                    maxupper=NULL,  # maximum for grid search of upper limit
                    minlower=NULL,  # minimum for grid search of lower limit
                    minupper=NULL,  # minimum for grid search of upper limit
                    maxlower=NULL){ # maximum for grid search of lower limit
    # estimation without constraint
    KMfit <- KMwithSE(tstar=tstar,data=data)
    KMest <- KMfit$res.tstar$Surv
    se.KM <- KMfit$res.tstar$se
    # min-max for searching for the upper and lower limit of CI
    if(missing(minupper)){
        minupper <- KMest + mytol
    }
    if(missing(maxupper)){
        maxupper <- 1- mytol
    }

    if(missing(minlower)){
        minlower <- mytol
    }
    if(missing(maxlower)){
        maxlower <- KMest - mytol
    }
    # upper
    upper <- NA
    try(upper <- stats::uniroot(f=function(x){LR.test.KM(tstar=tstar, pstar=x,  data=data)$p - alpha},
                                lower=minupper,
                                upper=maxupper, 
                                tol=mytol)$root, silent = TRUE)
    # lower
    lower <- NA
    try(lower <- stats::uniroot(f=function(x){LR.test.KM(tstar=tstar, pstar=x,  data=data)$p - alpha},
                                lower=minlower, 
                                upper=maxlower,
                                tol=mytol)$root, silent = TRUE)
    # return results
    CI <- c(lower,upper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}

#----------------------------------------------------------------------
### LRCIKM.R ends here
