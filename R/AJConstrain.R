AJConstrain <- function (tstar,
                         CIF1star,
                         data, mytol=10^(-5),
                         Trace=FALSE,
                         plotUNBUG = FALSE) {
    if(Trace) {
        print(paste("New call to AJConstrain, with CIF1star = ",CIF1star))
    }
    # {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    # }}}
    if(length(unique(time))!=length(time)){
        stop("This function does not handle ties to compute the constrained version of the Aalen-Johansen estimates.")
    }
    # {{{ useful quantities
    N <- length(time) 
    ## tab <- table(time, status)
    tab <- table(time,factor(status,levels=c(0,1,2)))
    d1i <- tab[, 2]      
    d2i <- tab[, 3]      
    ci <- tab[, 1]
    ti <- as.numeric(dimnames(tab)[[1]])
    k <- length(ti)
    ni <- c(N, N - cumsum(ci)[-k] - cumsum(d1i)[-k] - cumsum(d2i)[-k])
    # keep only for ti before tstar
    II <- ti <= tstar    
    nj <- ni[II]
    d1j <- d1i[II]
    d2j <- d2i[II]
    cj <- ci[II]
    kj <- sum(II)
    ## cbind(ti[II],nj,d1j,d2j,cj)
    # }}}
    # {{{ function to compute estimates given lamda
    CIFlambda <- function(lambda) {
        # initialize vectors
        a1 <- rep(0,kj) # cause-1-specific hazard
        a2 <- rep(0,kj) # cause-2-specific hazard
        a <-  rep(0,kj) # all-causes hazard
        Surv <- rep(1,kj) # survival : P(T > t_i)
        CIF1 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=1)
        CIF2 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=2)
        # loop over ti
        for (i in 1:kj) {
            if(i==1){
                if(d1j[i]>0){
                    # event at ti is event 1
                    a1[i] <- 1/(nj[i] + lambda*(1-0-CIF1star) )
                }else{
                    if(d2j[i]>0){
                        # event at ti is event 2
                        a2[i] <- 1/(nj[i] + lambda*(0-CIF1star))
                    }
                } # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
                a[i] <- a1[i] + a2[i]    
                Surv[i] <- 1*(1- a[i])
                CIF1[i] <- 0 + 1*a1[i]
                CIF2[i] <- 0 + 1*a2[i]                
            }else{
                if(d1j[i]>0){
                    # event at ti is event 1
                    a1[i] <-  1/(nj[i] + lambda*(1 - CIF2[i-1] - CIF1star) )
                }else{
                    if(d2j[i]>0){
                        # event at ti is event 2
                        a2[i] <- 1/(nj[i] + lambda*(CIF1[i-1] - CIF1star))
                    }
                }
                # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
                a[i] <- a1[i] + a2[i]    
                Surv[i] <- Surv[i-1]*(1-a[i])
                CIF1[i] <- CIF1[i-1] + Surv[i-1]*a1[i]
                CIF2[i] <- CIF2[i-1] + Surv[i-1]*a2[i]                               
            }
        }
        # output
        out <- list(forRoot=CIF1[kj],
                    Surv=Surv,
                    a=a,
                    a1=a1,
                    a2=a2,
                    CIF1=CIF1,
                    CIF2=CIF2,
                    alla1a2Positive=all(c(a1>=0,a2>=0)),
                    alla1Positive=all(a1>=0),
                    alla2Positive=all(a2>=0),
                    allOK=all(c(a1>=0,
                                a2>=0,
                                CIF1<=1,
                                CIF1>=0,
                                CIF2<=1,
                                CIF2>=0,
                                CIF1 + CIF2<=1
                                ))
                    )
        out
    }
    # }}}
    # {{{ find lambda
    ## browser()
    # {{{ we find interval for lambda with positive values for all aki in which to use uniroot
    if(Trace){print(paste("lambda.min : start"))}
    thelambdamin <- stats::uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                   lower=min(-N*1000,-N*(10/mytol)), # How to choose that ????
                                   upper=0,
                                   tol=mytol)$root
    if(!CIFlambda(thelambdamin)$allOK){
        nwhilellop <- 0
        while(!CIFlambda(thelambdamin)$allOK & nwhilellop<30){
            if(Trace){ print(paste("While loop (1), step ",nwhilellop,": thelambdamin=",thelambdamin))}
            nwhilellop <- nwhilellop + 1
            thelambdamin <- stats::uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                           lower=thelambdamin,
                                           upper=0,
                                           ## extendInt = c("downX"),
                                           tol=mytol)$root
        }
    }
    if(Trace){
        print(paste("lambda.min : done"))
        print(paste("thelambdamin=",thelambdamin))
    }
    ## browser()
    if(sum(data$status==2 & data$time <= tstar)>=1){ # only if there are competing events before t
        thelambdamax <- stats::uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                       lower=0,
                                       ## upper=N*1000, # How to choose that ????
                                       upper=max(N*1000,N*(10/mytol)),
                                       ## extendInt = c("upX"),
                                       tol=mytol)$root
        if(Trace){
            print(paste("lambda.max : first done"))
        }
        if(!CIFlambda(thelambdamax)$allOK){
            nwhilellop2 <- 0
            while(!CIFlambda(thelambdamax)$allOK & nwhilellop2<30){
                if(Trace){
                    print(paste("While loop (2), step ",nwhilellop2,": thelambdamax=",thelambdamax))
                }
                nwhilellop2 <- nwhilellop2 + 1
                thelambdamax <- stats::uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                               lower=0,
                                               upper=thelambdamax,
                                               tol=mytol)$root
                if(Trace){
                    print(paste("while loop :",nwhilellop2))
                }
            }
        }
    }else{
        thelambdamax <- N*1000
    }
    if(Trace){
        print(paste("thelambdamax=",thelambdamax))
        print(paste("lambda.max : done"))
    }
    # }}}
    # {{{ to understand bugs    
    if(plotUNBUG){
        lengthxxx <- 1000
        xxx <- seq(from=thelambdamin,to=thelambdamax,length.out=lengthxxx)
        yyy <- sapply(xxx,function(x){o <- CIFlambda(x);o$forRoot - CIF1star})
        graphics::plot.default(xxx,yyy,ylim=c(-1,1))
        graphics::abline(h=0)
        graphics::abline(v=thelambdamin)
        graphics::abline(v=thelambdamax)
        ## graphics::abline(v=(-N/(1-CIF1star)),col="red")
        ## graphics::abline(v=(N/CIF1star),col="red")
        #    
        yyy[1]
        yyy[lengthxxx]
        CIFlambda(xxx[1])$alla1a2Positive
        CIFlambda(xxx[lengthxxx])$alla1a2Positive
        CIFlambda(xxx[1] + 5)$alla1a2Positive
        CIFlambda(xxx[1] + 100/N)$forRoot - CIF1star
        CIFlambda(xxx[lengthxxx])$alla1a2Positive
        #
        CIFlambda(thelambdamin)$alla1a2Positive
        CIFlambda(thelambdamax)$alla1a2Positive
        CIFlambda(thelambdamax)$forRoot
        CIFlambda(thelambdamax)$forRoot - CIF1star
        CIFlambda(thelambdamin)$forRoot - CIF1star
        graphics::points(thelambdamin,CIFlambda(thelambdamin)$forRoot - CIF1star,col="red")
        graphics::points(thelambdamax,CIFlambda(thelambdamax)$forRoot - CIF1star,col="red")
        CIFlambda(thelambdamax)$allOK
        CIFlambda(thelambdamax)$CIF1
        CIFlambda(thelambdamax)$CIF2
        CIFlambda(thelambdamax)$a
    }
    # }}}
    if(Trace){print(paste("Int lambda.min , lambda.max : done"))}
    # {{{ find lambda within interval
    thelambda <- stats::uniroot(f=function(x){o <- CIFlambda(x);o$forRoot} - CIF1star,
                                lower=(thelambdamin + 1/N),
                                upper=(thelambdamax -1/N),tol=mytol)$root

    if(Trace){print("thelambda: done")}
    # {{{ to understand bugs
    if(plotUNBUG){
        graphics::abline(v=thelambda,col="red")
    }
    # }}}

    # }}}
    # }}}
    # {{{ compute everything with the right value for lambda
    outCIFlambda <- CIFlambda(thelambda)
    outCIFlambda0 <- CIFlambda(0)
    # }}}
    # {{{ output    
    out <- list(n=N,
                time = time[II],
                CIF1c = outCIFlambda$CIF1,
                CIF2c = outCIFlambda$CIF2,
                tstar=tstar,
                CIF1star=CIF1star,
                outCIFlambda=outCIFlambda,
                outCIFlambda0=outCIFlambda0,
                lambda=thelambda,
                lambdaSeachInt=c(thelambdamin,thelambdamax)
                )
    class(out) <- "constrainedAJ"
    # }}}
    if(Trace){
        print(paste("CIF1star = ",CIF1star,": DONE"))
    }
    out
}


