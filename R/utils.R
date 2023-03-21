### utils.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (12:04) 
## Version: 
## Last-Updated: Mar 21 2023 (15:13) 
##           By: Paul Blanche
##     Update #: 16
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

specdec <- function(x,k=2) {format(x=round(x=x,digits=k),nsmall=k)}


# {{{ lowerRisk and upperRisk: compute CI for risk using normal approx on the cloglog scale
lowerRisk <- function(est,se,alpha) {
    1-(1-est)^(exp(stats::qnorm(1-alpha/2)*se/((1-est)*log(1-est))))
}
upperRisk <- function(est,se,alpha) {
    1-(1-est)^(exp(-stats::qnorm(1-alpha/2)*se/((1-est)*log(1-est))))
}
# }}}

# {{{ ciRisk: Compute Wald-type confidence intervals for a risk ratio or a risk difference
# Compute Wald-type confidence intervals for a risk ratio or a risk difference
#
#  R1 pointwise estimate of the risk in group 1
#  R2 pointwise estimate of the risk in group 2 (the reference group)
#  seR1 se for risk in group 1
#  seR2 se for risk in group 2
#  alpha type-I error
#  RR.H0 risk ratio under the null hypothesis, for which a p-value is computed.
#  Diff.H0 risk difference under the null hypothesis, for which a p-value is computed.

ciRisk <- function(R1,R2,seR1,seR2,alpha=0.05,
                   RR.H0=1,
                   Diff.H0=0) {
    theq <- stats::qnorm(1-alpha/2)
    # {{{ Relative risk
    # CI computed with normal approximation of log(RR)
    RR <- R1/R2
    selogRR <- sqrt( (seR1^2)/(R1^2) + (seR2^2)/(R2^2)  )   # se of log(RR)
    lowerRR <- RR*exp(-theq*selogRR)
    upperRR <- RR*exp(theq*selogRR)
    zlogRR <- (log(RR) - log(RR.H0))/selogRR
    pRR <- 2*(1-stats::pnorm(abs(zlogRR)))
    allRR <- c(RR,selogRR,lowerRR,upperRR,pRR)
    # }}}    
    # {{{ Risk Difference
    # CI computed with normal approximation of the risk difference
    Rdiff <- R1-R2
    seRdiff <- sqrt(seR1^2 + seR2^2)
    lowerRdiff <- Rdiff - theq*seRdiff
    upperRdiff <- Rdiff + theq*seRdiff
    zRdiff <- (Rdiff-Diff.H0)/seRdiff
    pRdiff <- 2*(1-stats::pnorm(abs(zRdiff)))
    allRdiff <- c(Rdiff,seRdiff,lowerRdiff,upperRdiff,pRdiff)
    # }}}
    # {{{ CI of risks in each group
    # normal approximation of the risk
    NlowerR1 <- R1 - theq*seR1
    NupperR1 <- R1 + theq*seR1
    NlowerR2 <- R2 - theq*seR2
    NupperR2 <- R2 + theq*seR2
    # normal approximation of the cloglog(risk)
    lowerR1 <- lowerRisk(est=R1,se=seR1,alpha=alpha)
    upperR1 <- upperRisk(est=R1,se=seR1,alpha=alpha)
    lowerR2 <- lowerRisk(est=R2,se=seR2,alpha=alpha)
    upperR2 <- upperRisk(est=R2,se=seR2,alpha=alpha)
    allR1 <- rbind(c(R1,seR1,NlowerR1,NupperR1,NA),
                   c(R1,seR1,lowerR1,upperR1,NA))
    allR2 <- rbind(c(R2,seR2,NlowerR2,NupperR2,NA),
                   c(R2,seR2,lowerR2,upperR2,NA))
    rownames(allR1) <- rownames(allR2) <- c("plain","cloglog")                  
    # }}}
    # {{{ output
    names(allRdiff) <-  names(allRR) <- colnames(allR2) <- colnames(allR1) <- c("est.","se","lower","upper","pvalue")
    out <- list(R1=allR1,
                R2=allR2,
                RR=allRR,
                Rdiff=allRdiff)
    # }}}
    out
}
# }}}

#----------------------------------------------------------------------
### utils.R ends here
