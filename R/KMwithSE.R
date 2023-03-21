### KMwithSE.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jul  8 2022 (09:39) 
## Version: 
## Last-Updated: Mar 20 2023 (10:41) 
##           By: Paul Blanche
##     Update #: 28
#----------------------------------------------------------------------
## 
### Commentary:
##
##  Compute Kaplan-Meier estimator and associated Greenwood variance.
##  Very much like the code of the KMlambda function, which is strongly
##  inspired from the kmgw.calc function packge bpcp.
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

KMwithSE <- function (tstar,
                      data) 
{
    # {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    # }}}
    # {{{ From kmgw.calc function packge bpcp
    N <- length(time)
    tab <- table(time, status)
    dstat <- dimnames(tab)$status
    if (length(dstat) == 2) {
        di <- tab[, 2]      
        ci <- tab[, 1]
    }
    else if (length(dstat) == 1) {
        if (dstat == "1") {
            di <- tab[, 1]
            ci <- rep(0, length(tab))
        }
        else if (dstat == "0") {
            ci <- tab[, 1]
            di <- rep(0, length(tab))
        }
    }
    else stop("status should be 0 or 1")
    y <- as.numeric(dimnames(tab)[[1]])
    k <- length(y)
    ni <- c(N, N - cumsum(ci)[-k] - cumsum(di)[-k])
    names(ni) <- names(di)
    ni <- as.numeric(ni)
    di <- as.numeric(di)
    ci <- as.numeric(ci)
    KM <- cumprod((ni - di)/ni)
    gw <- KM^2 * cumsum(di/(ni * (ni - di)))
    gw[KM == 0] <- 0
    #--
    x <- list(time = y, # times 
              di = di,  # numbers of observed event at the times
              ni = ni,  # numbers of subjects at risk at the times
              KM = KM,  # KM estimates at the times
              gw = gw   # Greenwood estimate of the variance of KM estimates at the times
              )
    # }}}
    # {{{ Results at tstar
    II <- which(x$time <= tstar)
    var.gw.tstar <- gw[max(II)]
    Surv.gw.tstar <- KM[max(II)]
    # }}}
    # {{{ output
    out <- list(res.full=x,
                res.tstar=list(se=sqrt(var.gw.tstar),
                               Surv=Surv.gw.tstar)
                )
    # }}}    
    out
}   


#----------------------------------------------------------------------
### KMwithSE.R ends here
