### TwoGroupKMConst2.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (13:05) 
## Version: 
## Last-Updated: Mar 21 2023 (14:24) 
##           By: Paul Blanche
##     Update #: 144
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


TwoGroupKMConst2 <- function (tstar,
                              RR,
                              Rdiff,
                              SR,   # survival ratio
                              data,
                              mytol=10^(-5),gridpt=NULL){

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
    #--
    d1 <- data[data$group==1,]
    d0 <- data[data$group==0,]
    # }}}

    # {{{  Compute log-likelihood under constraint for either RR or Rdiff, for a specific t-year risk = pt for one group (here group 1)  
    ptToLogLikelihood0 <- function(pt){
        # start by group 1, then group 0
        # group 1
        res1 <- KMConstrain(tstar=tstar,
                            pstar=1-pt,
                            data=d1,
                            mytol=mytol)       
        #--
        if(!isMissRdiff){
            pt0 <- min(max(pt-Rdiff,mytol),1-mytol)
        }
        if(!isMissRR){
            pt0 <- min(pt/RR,1-mytol)
        }
        if(!isMissSR){
            pt0 <- max(1-(1-pt)/SR,mytol)
        }
        # group 0        
        res0 <- KMConstrain (tstar=tstar,
                             pstar=1-pt0,
                             data=d0,mytol=mytol)              
        ## Compute constrained log likelihood
        consLogLike0 <- MyLogLikeKM(a=res0$ac,
                                    ni=res0$nj,
                                    di=res0$dj)
        consLogLike1 <- MyLogLikeKM(a=res1$ac,
                                    ni=res1$nj,
                                    di=res1$dj)
        # output
        out <- c(consLogLike0=consLogLike0,
                 consLogLike1=consLogLike1,
                 lambda0=res0$lambda,
                 lambda1=res1$lambda)
        out
    }
    # Same but less in output: only likelihood under RR constraint, convenient for optimization.
    ptToLogLikelihood <- function(pt){
        x <- ptToLogLikelihood0(pt)
        -(x["consLogLike0"]+x["consLogLike1"])
    }   
    # }}}

    # {{{ Set upper and lower limit to search for pt  (i.e. constraint for group 1) in optim just below.
    # these lower and upper limits are just consequences that BOTH pt and pt0 should lie in (0,1) and that
    # we have the constraint on either RR, SR or Rdiff.
    if(!isMissRdiff){
        lowerOptim <- max(mytol,Rdiff+mytol)
        upperOptim <- min(1-mytol,1+Rdiff-mytol)       
    }
    if(!isMissRR){
        lowerOptim <- mytol
        upperOptim <- min(1-mytol,RR-mytol)
    }
    if(!isMissSR){
        lowerOptim <- max(mytol,1-SR+mytol)
        upperOptim <- 1-mytol
    }
    # }}}
    
    # {{{ optimization: search pt (i.e. constraint for group 1)
    res <- stats::optim(fn=ptToLogLikelihood,
                        method="Brent",
                        par=0.5,
                        lower=lowerOptim,
                        upper=upperOptim,
                        control=list(abstol=mytol,maxit=500)
                        )
    contRisk <- res$par
    ConsLogLike <- res$value
    EachConsLogLike <- ptToLogLikelihood0(res$par)
    # }}}
    # {{{ Compute the log-likelihood under constraint for a grid of values of pt.
    # Mainly for checking and/or pedagogical purpose.
    AllConsLogLike <- NULL
    if(!is.null(gridpt)){
        AllConsLogLike <- sapply(gridpt,ptToLogLikelihood0)
    }
    # }}}
    # {{{ Just for making the output of the function complete: redo the computation of
    # the constraint estimate outside optim.
    # group 1
    res1 <- KMConstrain(tstar=tstar,
                        pstar=1-contRisk,
                        data=d1)       
    # group 0
    if(!isMissRdiff){
        pt0 <- min(max(contRisk-Rdiff,mytol),1-mytol)
    }
    if(!isMissRR){
        pt0 <- min(contRisk/RR,1-mytol)
    }
    if(!isMissSR){
        pt0 <- max(1-(1-contRisk)/SR,mytol)

    }
    res0 <- KMConstrain (tstar=tstar,
                         pstar=1-pt0,
                         data=d0)
    # }}}
    # {{{ output
    out <- list(res0=res0,
                res1=res1,
                cRisk0=pt0,
                cRisk1=contRisk,
                urisk0=max(1-res0$S0), # unconstrained risk in group 0
                urisk1=max(1-res1$S0), # unconstrained risk in group 1
                ConsLogLike=-ConsLogLike,
                consLogLike0=EachConsLogLike["consLogLike0"],
                consLogLike1=EachConsLogLike["consLogLike1"],
                ResGrid=AllConsLogLike)
    # }}}
    out
}

#----------------------------------------------------------------------
### TwoGroupKMConst2.R ends here
