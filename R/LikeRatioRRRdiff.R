### LikeRatioRRRdiff.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 29 2018 (09:00) 
## Version: 
## Last-Updated: Aug 14 2023 (16:13) 
##           By: Paul Blanche
##     Update #: 94
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

LikeRatioRRRdiff <- function(data,
                             Rdiff,
                             RR,
                             tstar,
                             mytol=10^(-5),
                             whichAlgo=1
                             ){
    ## {{{ compute constrained CIF
    if((missing(Rdiff) & missing(RR)) | (!missing(Rdiff) & !missing(RR))){
        stop("Provide input for Either RR or Rdiff (but not both)")
    }else{
        if(!missing(RR)){
            if(whichAlgo==1){
                x <- TwoGroupAJConst(tstar=tstar,
                                     RR=RR,
                                     data=data,
                                     mytol=mytol)
            }
            if(whichAlgo==2){
                x <- TwoGroupAJConst2(tstar=tstar,
                                      RR=RR,
                                      data=data,
                                      mytol=mytol)
            }
            Rdiff <- NA
            unconstRdiff <- NA
            unconstRR <- x$unconstRR
        }else{
            if(whichAlgo==1){
                x <- TwoGroupAJConst(tstar=tstar,
                                     Rdiff=Rdiff,
                                     data=data,
                                     mytol=mytol)
            }
            if(whichAlgo==2){
                x <- TwoGroupAJConst2(tstar=tstar,
                                      Rdiff=Rdiff,
                                      data=data,
                                      mytol=mytol)
            }
            RR <- NA
            unconstRR <- NA
            unconstRdiff <- x$unconstRdiff
        }
    }
    ## }}}
    
    ## {{{ 
    resgroup1null <- CIFlambda(tstar=tstar,
                               lambda=0,
                               CIFstar=0.5,
                               data=data[data$group==1,],
                               Trace=FALSE)
    resgroup0null <- CIFlambda(tstar=tstar,
                               lambda=0,
                               CIFstar=0.5,
                               data=data[data$group==0,],
                               Trace=FALSE)
    ## }}}
    ## {{{ Compute constrained log likelihood
    consLogLike0 <- x$consLogLike0
    consLogLike1 <- x$consLogLike1
    ## }}}
    ## {{{ Compute unconstrained log likelihood
    LogLike0 <- MyLogLike(a1=resgroup0null$outCIFlambda$a1,
                          a2=resgroup0null$outCIFlambda$a2,
                          ssn=resgroup0null$n)
    LogLike1 <- MyLogLike(a1=resgroup1null$outCIFlambda$a1,
                          a2=resgroup1null$outCIFlambda$a2,
                          ssn=resgroup1null$n)
    ## }}}
    ## {{{ Likelihood ratio test and p-value
    LogLikRatio <- consLogLike0 + consLogLike1 - ( LogLike0 + LogLike1 )
    thepval <- 1-stats::pchisq(-2*(LogLikRatio),df=1)
    ## }}}
    out <- list(pval=thepval,
                RR=RR,
                Rdiff=Rdiff,
                unconstRR=unconstRR,
                unconstRdiff=unconstRdiff,
                LogLikRatio=LogLikRatio,
                consLogLike0=consLogLike0,
                consLogLike1=consLogLike1,
                LogLike0=LogLike0,
                LogLike1=LogLike1                
                )
    out
}

##----------------------------------------------------------------------
### LikeRatioRRRdiff.R ends here
