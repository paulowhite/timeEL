### KMlambda.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec 15 2021 (13:34) 
## Version: 
## Last-Updated: Mar 21 2023 (14:37) 
##           By: Paul Blanche
##     Update #: 11
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


KMlambda <- function (tstar, data, lambda) {
    # {{{ data management steps (just for me)
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
    I <- di > 0
    x <- list(time = y[I], # times (uncensored observations only)
              di = di[I],  # numbers of observed event at the times
              ni = ni[I],  # numbers of subjects at risk at the times
              KM = KM[I],  # KM estimates at the times
              gw = gw[I]   # Greenwood estimate of the variance of KM estimates at the times
              )
    # }}}
    # {{{ Very much inspired from kmConstrain of bpcp package
    II <- x$time <= tstar
    if (any(II)) {
        nj <- x$ni[II]
        dj <- x$di[II]
        pbar <- (nj + lambda - dj)/(nj + lambda)
        Sbar <- cumprod(pbar)
        # "hazard" estimate with contraint
        ac <-  dj/(nj + lambda)
        # add "hazard" estimate without contraint
        a0 <-  dj/nj
        S0 <- cumprod(1-a0)
        out <- list(time = x$time[II],
                    Sc = Sbar,
                    S0=S0,
                    lambda=lambda,
                    x=x,
                    tstar=tstar,
                    ac=ac,
                    a0=a0,
                    n=N,
                    nj=nj,
                    dj=dj
                    )
    }
    else {
        out <- list(time = tstar,
                    x=x,
                    tstar=tstar,
                    lambda=lambda,
                    S0=NA,
                    ac=NA,
                    a0=NA,
                    nj=NA,
                    dj=NA,
                    n=N
                    )
    }
    out
    # }}}    
}

#----------------------------------------------------------------------
### KMlambda.R ends here
