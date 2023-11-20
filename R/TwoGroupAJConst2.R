### TwoGroupAJNew.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 19 2020 (16:21) 
## Version: 
## Last-Updated: Aug 10 2023 (15:10) 
##           By: Paul Blanche
##     Update #: 268
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


TwoGroupAJConst2 <- function (tstar,
                              RR,
                              Rdiff,
                              data,
                              mytol=10^(-5),gridpt=NULL){

    # {{{ data management steps
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    if((missing(Rdiff) & missing(RR)) | (!missing(Rdiff) & !missing(RR))){
        stop("Provide input for Either RR or Rdiff (but not both)")
    }else{
        if(missing(RR)){
            isRdiff <- TRUE
        }else{
            isRdiff <- FALSE
        }
    }
    d1 <- data[data$group==1,]
    d0 <- data[data$group==0,]
    if(length(unique(d1$time))!=length(d1$time) | length(unique(d0$time))!=length(d0$time)){
        stop("This function does not handle ties to compute the constrained version of the Aalen-Johansen estimates.")
    }
    # }}}

    # {{{  Compute log-likelihood under constraint for either RR or Rdiff, for a specific t-year risk = pt for one group (here group 1)  
    ptToLogLikelihood0 <- function(pt){
        # start by group 1, then group 0
        # group 1
        res1 <- AJConstrain(tstar=tstar,
                            CIF1star=pt,
                            data=d1,
                            mytol=mytol,
                            Trace=FALSE)       
        # group 0
        pt0 <- ifelse(isRdiff,min(max(pt-Rdiff,mytol),1-mytol),min(pt/RR,1-mytol))
        ## print(paste0("pt0=",pt0))
        ## browser()
        res0 <- AJConstrain (tstar=tstar,
                             CIF1star=pt0,
                             data=d0,
                             mytol=mytol,
                             Trace=FALSE)              
        ## Compute constrained log likelihood
        consLogLike0 <- MyLogLike(a1=res0$outCIFlambda$a1,
                                  a2=res0$outCIFlambda$a2,
                                  ssn=res0$n)
        consLogLike1 <- MyLogLike(a1=res1$outCIFlambda$a1,
                                  a2=res1$outCIFlambda$a2,
                                  ssn=res1$n)
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
    # we have the constraint on either RR or Rdiff.
    if(isRdiff){
        lowerOptim <- max(mytol,Rdiff+mytol)
        upperOptim <- min(1-mytol,1+Rdiff-mytol)       
    }else{
        lowerOptim <- mytol
        upperOptim <- min(1-mytol,RR-mytol)
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
    res1 <- AJConstrain (tstar=tstar,
                         CIF1star=contRisk,
                         data=d1,
                         mytol=mytol,
                         Trace=FALSE)       
    # group 0
    pt0 <- ifelse(isRdiff,contRisk-Rdiff,contRisk/RR)
    res0 <- AJConstrain (tstar=tstar,
                         CIF1star=pt0,
                         data=d0,
                         mytol=mytol,
                         Trace=FALSE)
    # }}}
    # {{{ output
    out <- list(res0=res0,
                res1=res1,
                cRisk0=pt0,
                cRisk1=contRisk,
                ConsLogLike=-ConsLogLike,
                consLogLike0=EachConsLogLike["consLogLike0"],
                consLogLike1=EachConsLogLike["consLogLike1"],
                ResGrid=AllConsLogLike)
    # }}}
    out
}


    
#----------------------------------------------------------------------
### TwoGroupAJNew.R ends here
