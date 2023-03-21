

LRCIAJ <- function (tstar,
                      data,
                      alpha=0.05,
                      mytol=10^(-5),
                      maxupper=NULL, # maximum for grid search of upper limit
                      minlower=NULL, # minimum for grid search of lower limit
                      minupper=NULL, # minimum for grid search of upper limit
                      maxlower=NULL){ # maximum for grid search of lower limit
    # estimation without constraint
    res.AJwithSE <- AJwithSE(tstar=tstar,data=data,CompAllse=FALSE)
    CIF1star <- res.AJwithSE$risk.t
    # min-max for searching for the upper and lower limit of CI
    if(missing(minupper)){
        minupper <- CIF1star + mytol
    }
    if(missing(maxupper)){
        maxupper <- 1- mytol
    }

    if(missing(minlower)){
        minlower <- mytol
    }
    if(missing(maxlower)){
        maxlower <- CIF1star - mytol
    }
    # upper
    CIFupper <- NA
    try(CIFupper <- stats::uniroot(f=function(x){LR.test.AJ(tstar, x,  data)$p - alpha},
                                   lower=minupper,
                                   upper=maxupper, 
                                   tol=mytol)$root, silent = TRUE)
    # lower
    CIFlower <- NA
    try(CIFlower <- stats::uniroot(f=function(x){LR.test.AJ(tstar, x,  data)$p - alpha},
                                   lower=minlower, 
                                   upper=maxlower,
                                   tol=mytol)$root, silent = TRUE)
    # return results
    CI <- c(CIFlower,CIFupper)
    names(CI) <- paste(round(100*c(alpha/2,1-alpha/2),2),"%")
    return(CI)    
}
