### LikeRatioRRRdiffKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (15:39) 
## Version: 
## Last-Updated: Mar 21 2023 (15:06) 
##           By: Paul Blanche
##     Update #: 71
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


LikeRatioRRRdiffKM <- function(data,
                               Rdiff,
                               RR,
                               SR,
                               tstar,
                               mytol=10^(-5),
                               whichAlgo=2
                               ){
    if( sum(c(missing(Rdiff),missing(RR),missing(SR)))!=2 ){
        stop("Provide only one of RR, Rdiff or SR as input.")
    }
    unconstRdiff <- NA
    unconstSR <- NA
    unconstRR <- NA
    ## {{{ compute constrained KM
    if(!missing(RR)){
        x <- TwoGroupKMConst2(tstar=tstar,
                              RR=RR,
                              data=data,
                              mytol=mytol)
        unconstRR <- x$urisk1/x$urisk0
    }
    if(!missing(Rdiff)){
        x <- TwoGroupKMConst2(tstar=tstar,
                              Rdiff=Rdiff,
                              data=data,
                              mytol=mytol)
        unconstRdiff <- x$urisk1-x$urisk0
    }
    if(!missing(SR)){
        x <- TwoGroupKMConst2(tstar=tstar,
                              SR=SR,
                              data=data,
                              mytol=mytol)
        unconstSR <- (1-x$urisk1)/(1-x$urisk0)
    }
    #--
    if(missing(RR)) RR <- NA
    if(missing(Rdiff)) Rdiff <- NA
    if(missing(SR)) SR <- NA
    ## }}}
    ## {{{ Extract constrained log likelihood
    consLogLike0 <- x$consLogLike0
    consLogLike1 <- x$consLogLike1
    ## }}}
    ## {{{ Compute unconstrained log likelihood
    LogLike0 <- MyLogLikeKM(a=x$res0$a0,
                            ni=x$res0$nj,
                            di=x$res0$dj)
    LogLike1 <- MyLogLikeKM(a=x$res1$a0,
                            ni=x$res1$nj,
                            di=x$res1$dj)
    ## }}}
    ## {{{ Likelihood ratio test and p-value
    LogLikRatio <- consLogLike0 + consLogLike1 - ( LogLike0 + LogLike1 )
    thepval <- 1-stats::pchisq(-2*(LogLikRatio),df=1)
    ## }}}
    out <- list(pval=thepval,
                RR=RR,
                SR=SR,
                Rdiff=Rdiff,
                unconstRR=unconstRR,
                unconstSR=unconstSR,
                unconstRdiff=unconstRdiff,
                LogLikRatio=LogLikRatio,
                consLogLike0=consLogLike0,
                consLogLike1=consLogLike1,
                LogLike0=LogLike0,
                LogLike1=LogLike1                
                )
    out
}


#----------------------------------------------------------------------
### LikeRatioRRRdiffKM.R ends here
