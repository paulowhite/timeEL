### TwoGroupKMConst1.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec 15 2021 (13:16) 
## Version: 
## Last-Updated: Mar 21 2023 (14:27) 
##           By: Paul Blanche
##     Update #: 121
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


TwoGroupKMConst1 <- function (tstar,
                              RR,
                              Rdiff,
                              SR,   # survival ratio
                              data,
                              mytol=10^(-5)){
    # {{{ data management steps
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    isMissRdiff <- missing(Rdiff)
    isMissRR <- missing(RR)
    isMissSR <- missing(SR)
    if( sum(c(isMissRdiff,isMissRR,isMissSR))!=2 ){
        stop("Provide only one of RR, Rdiff or SR as input.")
    }
    if(isMissSR){stop("The method 1 works only for a given SR. This is because the equations of Thomas and Grunkemeier (1975) work in that case only. The equations that we should use instead for RR or Rdiff are too different, hence not implemeted.")}
    d1 <- data[data$group==1,]
    d0 <- data[data$group==0,]
    # }}}
    ## From Paper Thomas and Grunkemeier (1975),
    ## we see that it is similar to the one-sample case, see eq 7.5 and eq 2.3.   
    # {{{ 
    ToUniFunc <- function(lambda){
        ResKMlambda1 <- KMlambda(tstar=tstar,
                                 data=d1,
                                 lambda=lambda)
        ResKMlambda0 <- KMlambda(tstar=tstar,
                                 data=d0,
                                 lambda=-lambda) # -lambda, because same reasoning as in the competing risks setting with RR.
        out <- min(ResKMlambda1$Sc)/min(ResKMlambda0$Sc) - SR
        out
    }
    # Limits for searching lambda
    Temp1 <- KMlambda(tstar=tstar,
                      data=d1,
                      lambda=0)
    Temp0 <- KMlambda(tstar=tstar,
                      data=d0,
                      lambda=0)
    #--    
    LowerLambda <- max(-min(Temp1$nj-Temp1$dj)+mytol,-1000*Temp0$nj[1]+mytol)
    UpperLambda <- min( min(Temp0$nj-Temp0$dj)-mytol, 1000*Temp1$nj[1])  
    # Searching lambda
    ToUniFunc(LowerLambda)
    ToUniFunc(UpperLambda)
    lambda <- stats::uniroot(ToUniFunc, lower=LowerLambda,upper=UpperLambda,tol=mytol)$root
    # }}}    
    
    # {{{
    ResKMlambda1 <- KMlambda(tstar=tstar,
                             data=d1,
                             lambda=lambda)
    ResKMlambda0 <- KMlambda(tstar=tstar,
                             data=d0,
                             lambda=-lambda)
    # }}}

    # {{{ Compute constrained log likelihood
    consLogLike0 <- MyLogLikeKM(a=ResKMlambda0$ac,
                                ni=ResKMlambda0$nj,
                                di=ResKMlambda0$dj)
    consLogLike1 <- MyLogLikeKM(a=ResKMlambda1$ac,
                                ni=ResKMlambda1$nj,
                                di=ResKMlambda1$dj)
    # }}}
   
    # {{{ output
    out <- list(resKM0=ResKMlambda0,
                resKM1=ResKMlambda1,
                cRisk0=max(1-ResKMlambda0$Sc),
                cRisk1=max(1-ResKMlambda1$Sc),
                urisk0=max(1-ResKMlambda0$S0), # unconstrained risk in group 0
                urisk1=max(1-ResKMlambda1$S0), # unconstrained risk in group 1)
                ConsLogLike=consLogLike0+consLogLike1,
                consLogLike0=consLogLike0,
                consLogLike1=consLogLike1)
    # }}}
    out
}


#----------------------------------------------------------------------
### TwoGroupKMConst1.R ends here
