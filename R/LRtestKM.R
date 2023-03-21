### Likelihood-hypo-test-KM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 30 2021 (15:10) 
## Version: 
## Last-Updated: Mar 21 2023 (15:06) 
##           By: Paul Blanche
##     Update #: 69
#----------------------------------------------------------------------
## 
### Commentary: 
##  Function to compute likelihood ratio, corresponding to comparison
##  between with and without constraint for Kaplan-Meier estimates.
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


LR.KM <- function(a.1,
                  a.2,
                  ni,
                  di
                  ){
    # loglikelihood corresponding to KM1
    ll1 <- MyLogLikeKM(a.1,ni,di)
    # loglikelihood corresponding to KM2
    ll2 <- MyLogLikeKM(a.2,ni,di)
    # difference in log-like
    LR <- ll2 - ll1
    LR
}


LR.test.KM <- function (tstar, pstar,  data, mytol=10^(-5)){
    # fit with constaint
    KM.cont.fit <- KMConstrain(tstar, pstar,  data, mytol)
    # log-like comparison
    myLR01 <- LR.KM(a.1 = KM.cont.fit$a0,
                    a.2 = KM.cont.fit$ac,
                    ni = KM.cont.fit$nj,
                    di = KM.cont.fit$dj
                    )
    # p-value from Chi-square distribution
    thepval <- 1-stats::pchisq(-2*myLR01,df=1)
    # output
    return(list(p=thepval,
                Surv.H0=pstar,
                Surv.Est=min(KM.cont.fit$S0),
                n=KM.cont.fit$n,
                TestStat=-2*myLR01))
}
#----------------------------------------------------------------------
### Likelihood-hypo-test-KM.R ends here
