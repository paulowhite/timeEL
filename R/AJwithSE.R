### AJwithSE.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec 27 2021 (11:30) 
## Version: 
## Last-Updated: Apr  4 2023 (20:54) 
##           By: Paul Blanche
##     Update #: 100
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


AJwithSE <- function (tstar,
                      data,CompAllse=FALSE) 
{
    # {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    # }}}
    # {{{ useful quantities
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
    # compute AJ at all times betore tstar
    a1 <- d1j/nj
    a2 <- d2j/nj
    a <- a1 + a2
    Surv <- c(1,cumprod(1- a))
    names(Surv)[1] <- 0
    CIF1 <- c(0,cumsum(Surv[-(kj+1)]*a1))
    CIF2 <- c(0,cumsum(Surv[-(kj+1)]*a2))
    risk.t <- CIF1[length(CIF1)]
    names(risk.t) <- NULL
    names(CIF1) <- names(CIF2) <- names(Surv)
    which((CIF1 + CIF2 + Surv)!=1)
    unique(CIF1 + CIF2 + Surv)    
    #--
    if(CompAllse){
        #--- compute se at all times s such that s<=t 
        se <- rep(NA,kj)
        for(k in 1:kj){
            line1 <- (((CIF1[k+1]- CIF1[2:(k+1)] )^2)/(nj[1:k] - (d1j[1:k] + d2j[1:k])))*(a1[1:k] + a2[1:k])
            line2 <- ((Surv[1:k]^2)/(nj[1:k]^3))*(d1j[1:k])*( nj[1:k] - d1j[1:k] - 2*(nj[1:k] - (d1j[1:k] + d2j[1:k]))*( CIF1[k+1]- CIF1[2:(k+1)] )/Surv[2:(k+1)])
            se[k] <- sqrt(sum(line1 + line2))        
        }
        se <- c(0,se)
    }else{
        #--- compute se at t only
        line1 <- (((max(CIF1)- CIF1[-1] )^2)/(nj - (d1j + d2j)))*(a1 + a2)
        line2 <- ((Surv[-(kj+1)]^2)/(nj^3))*(d1j)*( nj - d1j - 2*(nj - (d1j + d2j))*(max(CIF1)- CIF1[-1] )/Surv[-1]  )
        # {{{ To handle cases where
        # (nj - (d1j + d2j))
        # goes to zero
        line1 <- line1[!is.nan(line1)]
        line2 <- line2[!is.nan(line2)]
        # }}}
        
        se <- c(0,rep(NA,kj-1),sqrt(sum(line1 + line2)))
        #---
    }
    # output    
    out <- list(Surv=Surv,
                CIF=CIF1,
                se=se,
                n.risk=c(N,nj),
                times=as.numeric(names(Surv)),
                n.event.1=c(0,d1j),
                n.event.2=c(0,d2j),
                #---                
                risk.t=risk.t,
                se.t=se[length(se)]
                )
    out
}   

#----------------------------------------------------------------------
### AJwithSE.R ends here
