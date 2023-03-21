# {{{ LR.AJ :  V2
# Function to compute the log likelihood ratio (LR)

# input :
# - a1.1 : vector of cause specific hazards of cause 1 for CIF1
# - a2.1 : vector of cause specific hazards of cause 2 for CIF1
# - a1.2 : vector of cause specific hazards of cause 1 for CIF2
# - a2.2 : vector of cause specific hazards of cause 2 for CIF2
# ssn : total Sample Size n of the dataset (i.e. number at risk at time 0)

# output :
# - LR = l(\lambda_2) - l(\lambda_1) :  (Log of the) Likelihood Ratio

LR.AJ <- function(a1.1,
                  a2.1,
                  a1.2,
                  a2.2,
                  ssn
                  ){
    # indices a1 not equal to zeros (same with and without constraint)
    ssn <- ssn-1
    i1 <- which(a1.1 !=0)
    i2 <- which(a2.1 !=0)
    # number of observerved times before the timepoint of ineterest
    nbt <- length(a1.1)
    # loglikelihood corresponding to AJ1
    ll1 <- sum(log(a1.1[i1])) + sum(log(a2.1[i2])) + sum((ssn:(ssn-nbt+1))*log(1 - a1.1 - a2.1))
    # loglikelihood corresponding to AJ2
    ll2 <- sum(log(a1.2[i1])) + sum(log(a2.2[i2])) + sum((ssn:(ssn-nbt+1))*log(1 - a1.2 - a2.2))
    LR <- ll2 - ll1
    LR
}
# }}}


# {{{ LR.test.AJ
# input :
#
# - tstar : time point t of interest
# - CIF1star : value of CIF at t under H0
# - data
#
# output :
# - p : pvalue of LRT
# - CIFest : CIF Aalen-Johansen estimate (nonparametric without any constrainst
# - CIF0 = CIF1star : value of CIF at t under H0
# - n : sample size of the data

LR.test.AJ <- function (tstar, CIF1star,  data){
    AJ.cont.fit <- AJConstrain(tstar, CIF1star,  data)
    myLR01 <- LR.AJ(a1.1 = AJ.cont.fit$outCIFlambda0$a1,
                    a2.1 = AJ.cont.fit$outCIFlambda0$a2,
                    a1.2 = AJ.cont.fit$outCIFlambda$a1,
                    a2.2 = AJ.cont.fit$outCIFlambda$a2,
                    ssn = AJ.cont.fit$n)
    thepval <- 1-stats::pchisq(-2*myLR01,df=1)
    return(list(p=thepval,
                CIF0=CIF1star,
                CIFest=max(AJ.cont.fit$outCIFlambda0$CIF1),
                n=AJ.cont.fit$n,
                TestStat=-2*myLR01))
}
# }}}
