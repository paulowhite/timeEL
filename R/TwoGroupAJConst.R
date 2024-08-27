## This current version assumes there are no ties (times, within treatment subgroups)
                                        
TwoGroupAJConst <- function (tstar,
                             RR,
                             Rdiff,
                             data,
                             mytol=10^(-5),
                             Trace=FALSE,
                             ToDebug=FALSE,
                             ShowPlot=FALSE) 
{
    # {{{ data management steps
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    if((missing(Rdiff) & missing(RR)) | (!missing(Rdiff) & !missing(RR))){
        stop("Provide input for Either RR or Rdiff (but not both)")
    }else{
        if(missing(RR)){
            mycst <- Rdiff
            isRdiff <- TRUE
        }else{
            mycst <- 0
            isRdiff <- FALSE
        }
    }
    d1 <- data[data$group==1,]
    d0 <- data[data$group==0,]
    if(length(unique(d1$time))!=length(d1$time) | length(unique(d0$time))!=length(d0$time)){
        stop("This function does not handle ties to compute the constrained version of the Aalen-Johansen estimates.")
    }
    # }}}
    # {{{ unconstraint estimates
    fitAJ0 <- AJwithSE(tstar=tstar,data=d0,CompAllse=FALSE)
    fitAJ1 <- AJwithSE(tstar=tstar,data=d1,CompAllse=FALSE)   
    CIF1gr1 <- fitAJ1$risk.t
    CIF1gr0 <- fitAJ0$risk.t
    # }}}
    
    # {{{ subfunction
    # functon that as pi1 (i.e. constraint fot the risk at t) as input
    # and the corresponding lambda and pi1 as output
    Pi1ToResGroup1 <- function(pi1){
        resgroup1 <- AJConstrain (tstar=tstar,
                                  CIF1star=pi1,
                                  data=d1,
                                  mytol=mytol,
                                  Trace=Trace)
        out <- list(lambda=resgroup1$lambda,
                    pi1=pi1,
                    outCIFlambda=resgroup1$outCIFlambda,
                    n=resgroup1$n)
        out
    }
    # functon that as pi1 (or equivalently CIFstar,
    # i.e. constraint fot the risk at t) and lambda as input
    # and the corresponding  estimates of all CSH as output
    Group1ResToGroup0Res <- function(lambda,CIFstar){
        resgroup0 <- CIFlambda(tstar=tstar,
                               lambda=-lambda,
                               CIFstar=CIFstar,
                               data=d0,
                               Trace=Trace)
        resgroup0
    }
    FindPi1 <- function(data){
        mylowerpi1 <- mytol
        myupperpi1 <- 1-mytol
        #--
        bisfpi1 <- function(pi1){
            x1 <- Pi1ToResGroup1(pi1=pi1);
            x2 <- Group1ResToGroup0Res(lambda=ifelse(isRdiff,x1$lambda,RR*x1$lambda),
                                       CIFstar=ifelse(isRdiff,pi1-Rdiff,pi1/RR));
            ## mydiff <- ifelse(isRdiff,pi1-max(x2$CIF1c)-Rdiff,pi1/max(x2$CIF1c)-1);
            ifelse(isRdiff,pi1-Rdiff-max(x2$CIF1c),pi1-RR*max(x2$CIF1c))
        }
        # {{{ 
        # {{{ To find search interval (which does not lead to negative hazard)
        areAllCSHPos <- function(pi1){
            x1 <- Pi1ToResGroup1(pi1=pi1);
            x2 <- Group1ResToGroup0Res(lambda=ifelse(isRdiff,x1$lambda,RR*x1$lambda),
                                       CIFstar=ifelse(isRdiff,pi1-Rdiff,pi1/RR));
            x2$outCIFlambda$allOK        
        }
        if(ToDebug){print(paste("pi1t.min : start"))}
        thepi1tmin <- mylowerpi1
        try(thepi1tmin <- stats::uniroot(f=function(l){as.numeric(areAllCSHPos(l)) - 0.5},
                                         lower=mylowerpi1, 
                                         upper=CIF1gr1,
                                         tol=mytol)$root, silent = TRUE)
        if(ToDebug){
            print(paste("pi1t.min : done"))
            print(paste("pi1t.min=",round(thepi1tmin*100,3),"%"))
        }
        #-
        thepi1tmax <- myupperpi1
        try(thepi1tmax <- stats::uniroot(f=function(l){as.numeric(areAllCSHPos(l)) - 0.5},
                                         lower=CIF1gr1,
                                         upper=myupperpi1,
                                         tol=mytol)$root, silent = TRUE)
        if(ToDebug){
            print(paste("thepi1tmax : done"))
            print(paste("thepi1tmax=",round(thepi1tmax*100,3),"%"))
        }

        # }}}        
        if(ToDebug){
            xxx <- seq(from=mylowerpi1,
                       to=myupperpi1,
                       length.out=100)
            zzz <- sapply(xxx,bisfpi1)
            # {{{ to NOT change the user's options
            oldpar <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(oldpar)) 
            # }}}
            graphics::par(mfrow=c(1,3))
            graphics::plot.default(xxx,zzz,xlab="pi1", ylab="f for which we need the root")
            graphics::abline(v=CIF1gr1,col="blue",lty=2,lwd=2)
            #
            graphics::abline(v=thepi1tmax,lty=2,col="green", lwd=2)
            graphics::abline(v=thepi1tmin,lty=2,col="green", lwd=2)
            #
            graphics::abline(h=0,lty=2,col="red")
            xxxRoot <- xxx[xxx <= thepi1tmax & xxx >= thepi1tmin ][which.min(abs(zzz[xxx <= thepi1tmax & xxx >= thepi1tmin ]))]
            graphics::abline(v=xxxRoot,lty=2,lwd=2)
            ## graphics::text(x=xxx[which.min(abs(zzz))],
            ## y=CIF1gr1,
            ## labels=paste0("pi1=",round(xxx[which.min(abs(zzz))],3)),
            ## pos=1)
            graphics::legend("bottomleft",lty=2,lwd=2,col=c("blue","black","green"),
                             legend=c("Unconstraint","constraint","limits of search interval (where all Hazards > 0)"),
                             bty="n")
        }
        # }}}
        stats::uniroot(f=bisfpi1,
                       lower=thepi1tmin,
                       upper=thepi1tmax,
                       ## extendInt = c("downX"),
                       maxiter=5000,
                       ## tol=mytol)$root
                       tol=mytol)
    }
    # }}}
    # {{{ Find Pi1 (which uniquely defines everything)
    ResfindPi1 <- FindPi1(data)
    pi1 <- ResfindPi1$root
    # }}}
    # {{{ output
    resGroup1 <- Pi1ToResGroup1(pi1=pi1)
    resGroup0 <- Group1ResToGroup0Res(lambda=ifelse(isRdiff,resGroup1$lambda,RR*resGroup1$lambda),
                                      CIFstar=ifelse(isRdiff,pi1-Rdiff,pi1/RR)
                                      )
    lambda <- resGroup1$lambda
    if(ShowPlot){
        if(!ToDebug){
            # {{{ to NOT change the user's options
            oldpar <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(oldpar)) 
            # }}}            
            graphics::par(mfrow=c(1,2))
        }
        ResTemps1 <- AJConstrain (tstar=tstar,
                                  CIF1star=pi1,
                                  data=d1,
                                  mytol=mytol,
                                  Trace=FALSE)
        plot.constrainedAJ(ResTemps1)
        graphics::abline(h=CIF1gr1,col="black",lty=2,lwd=2)
        graphics::text(x=3*tstar/4,
                       y=CIF1gr1,col="black",
                       labels=paste0("No constraint=",round(CIF1gr1,3)),
                       pos=3)
        graphics::text(y=pi1,
                       x=tstar/4,
                       labels=paste0("Constrained=",round(pi1,3)),
                       pos=1)
        #--
        graphics::title(paste0("With constraint ",
                               ifelse(isRdiff,"Rdiff=","RR="),
                               ifelse(isRdiff,Rdiff,RR),"\n without constraint ",
                               ifelse(isRdiff,"Rdiff=","RR="),
                               round(ifelse(isRdiff,CIF1gr1-CIF1gr0,CIF1gr1/CIF1gr0),1),
                               "\n (group 1: exposed)"
                               ))
        plot.AJlambda(resGroup0)
        graphics::text(y=ifelse(isRdiff,pi1-Rdiff,pi1/RR),
                       x=3*tstar/4,
                       labels=paste0("Constrained=",round(ifelse(isRdiff,pi1-Rdiff,pi1/RR),3)),
                       pos=1)
        graphics::lines(fitAJ0$times,fitAJ0$CIF,type="s",col="red",
                        lty=2,lwd=2)
        d0comp <- d0
        d0comp$status[d0$status==1] <- 2
        d0comp$status[d0$status==2] <- 1       
        fitAJ0comp <- AJwithSE(tstar=tstar,data=d0comp,CompAllse=FALSE)
        graphics::lines(fitAJ0comp$times,fitAJ0comp$CIF,type="s",col="blue",
                        lty=2,lwd=2)       
        graphics::abline(h=CIF1gr0,col="black",lty=2,lwd=2)
        graphics::text(x=tstar/4,
                       y=CIF1gr0,col="black",
                       labels=paste0("No constraint=",round(CIF1gr0,3)),
                       pos=3)
        graphics::title(paste0("With constraint ",
                               ifelse(isRdiff,"Rdiff=","RR="),
                               ifelse(isRdiff,Rdiff,RR),"\n without constraint ",
                               ifelse(isRdiff,"Rdiff=","RR="),
                               round(ifelse(isRdiff,CIF1gr1-CIF1gr0,CIF1gr1/CIF1gr0),1),
                               "\n (group 0: non exposed)"
                               ))
        graphics::legend("left",
                         col=c("red","blue"),
                         lty=2,
                         lwd=2,
                         bty="n",
                         legend=c("AJ cause=1",
                                  "AJ cause=2")
                         )
    }

    # {{{ Constraint Log likelihood
    consLogLike0 <- MyLogLike(a1=resGroup0$outCIFlambda$a1,
                              a2=resGroup0$outCIFlambda$a2,
                              ssn=resGroup0$n)
    consLogLike1 <- MyLogLike(a1=resGroup1$outCIFlambda$a1,
                              a2=resGroup1$outCIFlambda$a2,
                              ssn=resGroup1$n)
    # }}}    
    # {{{ output
    out <- list(pi1=pi1,
                pi0=max(resGroup0$CIF1c),
                lambda=lambda,
                resGroup0=resGroup0,
                resGroup1=resGroup1,
                RR=ifelse(isRdiff,NA,RR),
                Rdiff=ifelse(isRdiff,Rdiff,NA),
                unconstRR=CIF1gr1/CIF1gr0,
                unconstRdiff=CIF1gr1-CIF1gr0,
                AbsRisk1=CIF1gr1,
                AbsRisk0=CIF1gr0,
                consLogLike0=consLogLike0,
                consLogLike1=consLogLike1
                )
    # }}}    

    # }}}
    out
}

