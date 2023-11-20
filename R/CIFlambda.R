### CIFlambda.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 23 2018 (09:03) 
## Version: 
## Last-Updated: Aug 14 2023 (16:10) 
##           By: Paul Blanche
##     Update #: 64
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

CIFlambda <- function (tstar,
                       lambda,
                       CIFstar,
                       data,
                       Trace=FALSE) 
{
    if(Trace){
        print(paste("New call to CIFlambda, with lambda = ",lambda))
    }
    ## {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    ## }}}
    if(length(unique(time))!=length(time)){
        stop("This function does not handle ties to compute the constrained version of the Aalen-Johansen estimates.")
    }
    ## {{{ useful quantities
    N <- length(time) 
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
    ## }}}
    ## {{{ Compute estimates given lamda

    ## {{{ initialize vectors
    a1 <- rep(0,kj) # cause-1-specific hazard
    a2 <- rep(0,kj) # cause-2-specific hazard
    a <-  rep(0,kj) # all-causes hazard
    Surv <- rep(1,kj) # survival : P(T > t_i)
    CIF1 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=1)
    CIF2 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=2)
    ## }}}
    ## {{{ loop over ti       
    for (i in 1:kj) {
        if(i==1){
            if(d1j[i]>0){
                # event at ti is event 1
                a1[i] <- 1/(nj[i] + lambda*(1-0-CIFstar) )
            }else{
                if(d2j[i]>0){
                    # event at ti is event 2
                    a2[i] <- 1/(nj[i] + lambda*(0-CIFstar))
                }
            } # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
            a[i] <- a1[i] + a2[i]    
            Surv[i] <- 1*(1- a[i])
            CIF1[i] <- 0 + 1*a1[i]
            CIF2[i] <- 0 + 1*a2[i]                
        }else{
            if(d1j[i]>0){
                # event at ti is event 1
                a1[i] <-  1/(nj[i] + lambda*(1 - CIF2[i-1] - CIFstar) )
            }else{
                if(d2j[i]>0){
                    # event at ti is event 2
                    a2[i] <- 1/(nj[i] + lambda*(CIF1[i-1] - CIFstar))
                }
            }
            # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
            a[i] <- a1[i] + a2[i]    
            Surv[i] <- Surv[i-1]*(1-a[i])
            CIF1[i] <- CIF1[i-1] + Surv[i-1]*a1[i]
            CIF2[i] <- CIF2[i-1] + Surv[i-1]*a2[i]                               
        }
    }
    ## }}}
    ## {{{ output
    out1 <- list(forRoot=CIF1[kj],
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
                             CIF1 + CIF2<=1+0.5*(1e-5)
                             )))
    out <- list(n=N,
                time = time[II],
                CIF1c = out1$CIF1,
                CIF2c = out1$CIF2,
                tstar=tstar,
                CIFstar=CIFstar,
                outCIFlambda=out1,
                lambda=lambda,
                lambdaSeachInt=c(NA,NA)
                )
    class(out) <- "AJlambda"
    ## }}}
    out
}


##----------------------------------------------------------------------
### CIFlambda.R ends here
