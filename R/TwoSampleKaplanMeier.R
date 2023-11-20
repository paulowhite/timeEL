### TwoSampleKaplanMeier.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (15:05) 
## Version: 
## Last-Updated: Aug 10 2023 (15:25) 
##           By: Paul Blanche
##     Update #: 265
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @description
##' Computes a risk difference, risk ratio or survival ratio with right-censored data, together with a confidence interval and
##' a p-value (to test for a difference between two groups). Pointwise estimates are computed via the Kaplan-Meier estimator. Computation of confidence intervals
##' and p-values are based on either Empirical Likelihood (EL) inference or Wald-type inference. Both are non-parametric approaches, which are asymptotically equivalent.
##' See Thomas & Grunkemeier (1975) for details about the Empirical Likelihood method. For the Wald-type approach, the asymptotic normal approximation is used on
##' the log scale for the risk ratio or survival ratio. No transformation is used for the risk or survival difference.
##'
##'
##' 
##' @title Risk difference and risk ratio using the Kaplan-Meier method
##' 
##' @param time vector of times (possibly censored)
##' @param status vector of usual survival status indicators (0 for censored observations, 1 otherwise)
##' @param group vector of binary group indicator. The reference group should be coded 0, the other 1. 
##' @param t the time point of interest (e.g. 1 to compute 1-year risk ratio)
##' @param SR.H0 the survival ratio under the null hypothesis, to compute a p-value. Default is 1.
##' @param RR.H0 the risk ratio under the null hypothesis, to compute a p-value. Default is 1.
##' @param Diff.H0 the risk difference under the null hypothesis, to compute a p-value. Default is 0. 
##' @param level confidence level for the confidence intervals. Default is 0.95. 
##' @param contr list of control parameters. tol=tolerance for numerical computation, default is 1e-5. method="EL", "Wald" or "both"  indicates wether 95% CI and the p-value should be computed based on Empirical Likelihood inference, Wald-type inference or both. algo=2 is currently the only option that is implemented.
##'
##' @return an object of class 'TwoSampleKaplanMeier'
##' 
##' @author Paul Blanche
##'
##' @references
##' Thomas & Grunkemeier (1975). Confidence interval estimation of survival probabilities for censored data. Journal of the American Statistical Association, 70(352), 865-871.
##' 
##' @examples
##' # This example reproduces some results presented in Table 4 of Thomas and Grunkemeier (1975)
##' Res2SKM95 <- TwoSampleKaplanMeier(time=Freireich$time,
##'                                   status=Freireich$status,
##'                                   group=Freireich$group,
##'                                   t=10)
##' Res2SKM95
##'
##'
##' @export
TwoSampleKaplanMeier <- function(time,
                                 status,
                                 group,
                                 t,
                                 SR.H0=1,
                                 RR.H0=1,
                                 Diff.H0=0,
                                 level=0.95,
                                 contr=list(tol=1e-5,algo=2,k=3,Trace=FALSE,method="both")){   

    # {{{ About input parameters
    # to handle the fact that we might want to give only some of the arguments in contr
    if(!missing(contr)){
        if(is.null(contr$tol))  contr$tol <- 1e-5       
        if(is.null(contr$algo)) contr$algo <- 2
        if(is.null(contr$k)) contr$k <- 3
        if(is.null(contr$Trace)) contr$Trace <- FALSE
        if(is.null(contr$method)) contr$method <- "both"
    }
    # basic checks
    if(!(contr$method %in% c("both","Wald","EL"))){
        stop("Choose contr$method='both' or 'Wald' or 'EL'.")
    }
    if(contr$algo!=2){
        stop("Choose algo=2. Currently this is the only option that is implemented.")
    }
    if(!("numeric" %in% class(time))){
        stop("time should be numeric.")
    }
    if(!all(status %in% c(0,1))){
        stop("status should contain values which are either 0 or 1.")
    }
    if(!all(group %in% c(0,1))){
        stop("group should contain values which are either 0 or 1.")
    }
    if(t <=0 ){
        stop("The time point of interest, t, should be positive.")
    }
    if(RR.H0<=0){
        stop("RR.H0 should be positive.")
    }
    if(SR.H0<=0){
        stop("SR.H0 should be positive.")
    }
    if(!(Diff.H0 > -1 & Diff.H0 < 1)){
        stop("Diff.H0 should be within ]-1,1[.")
    }
    if( !(level > 0 & level < 1)   ){
        stop("level should be within ]0,1[.")
    }
    # }}}


    # {{{ data 
    d <- data.frame(time=time,status=status,group=group)
    d0 <- d[d$group==0,]
    d1 <- d[d$group==1,]
    if(min(d0$time[d0$status==1])>t | min(d1$time[d1$status==1])>t){
        stop("There is no event observed in the data before time t in at least one of the two groups.")
    }
    # }}}

    # {{{ Wald inference
    KMfit0 <- KMwithSE(tstar=t,data=d0)
    KMfit1 <- KMwithSE(tstar=t,data=d1)
    #--
    risk1 <- 1-KMfit1$res.tstar$Surv
    risk0 <- 1-KMfit0$res.tstar$Surv
    se.risk1 <- KMfit1$res.tstar$se
    se.risk0 <- KMfit0$res.tstar$se
    surv1 <- 1-risk1
    surv0 <- 1-risk0   
    #--
    res.Wald <- ciRisk(R1=risk1,
                       R2=risk0,
                       seR1=se.risk1,
                       seR2=se.risk0,
                       RR.H0=RR.H0,
                       Diff.H0=Diff.H0,
                       alpha=1-level)
    res.Wald.SR <- ciRisk(R1=surv1,
                          R2=surv0,
                          seR1=se.risk1,
                          seR2=se.risk0,
                          RR.H0=RR.H0,
                          Diff.H0=Diff.H0,
                          alpha=1-level)
    # }}}
   
    if(contr$method=="both" | contr$method=="EL"){
        
        # {{{ Additional Wald inference, just to fine tune initial values for EL computation
        k <- contr$k
        outWald1 <- ciRisk(R1=risk1,
                           R2=risk0,
                           seR1=se.risk1,
                           seR2=se.risk0,
                           RR.H0=RR.H0,
                           Diff.H0=Diff.H0,
                           alpha=(1-level)/(k**1.3))
        #--
        outWald2 <- ciRisk(R1=risk1,
                           R2=risk0,
                           seR1=se.risk1,
                           seR2=se.risk0,
                           RR.H0=RR.H0,
                           Diff.H0=Diff.H0,
                           alpha=(1-level)*k)
        #--
        outWald1S <- ciRisk(R1=surv1,
                            R2=surv0,
                            seR1=se.risk1,
                            seR2=se.risk0,
                            RR.H0=RR.H0,
                            Diff.H0=Diff.H0,
                            alpha=(1-level)/(k**1.3))
        #--
        outWald2S <- ciRisk(R1=surv1,
                            R2=surv0,
                            seR1=se.risk1,
                            seR2=se.risk0,
                            RR.H0=RR.H0,
                            Diff.H0=Diff.H0,
                            alpha=(1-level)*k)
        # }}}
        # {{{ EL CI inference

        # {{{ Computation with first search interval
        outEmpLikeRR <- EmpLikRRCIKM(data=d,
                                     tstar=t,
                                     mytol=contr$tol,
                                     alpha=1-level,
                                     maxupper=outWald1$RR["upper"],
                                     minlower=outWald1$RR["lower"],
                                     minupper=outWald2$RR["upper"],
                                     maxlower=outWald2$RR["lower"],                        
                                     whichAlgo=contr$algo)
        #--
        outEmpLikeRdiff <- EmpLikRdiffCIKM(data=d,
                                           tstar=t,
                                           mytol=contr$tol,
                                           alpha=1-level,
                                           maxupper=outWald1$Rdiff["upper"],
                                           minlower=outWald1$Rdiff["lower"],
                                           minupper=outWald2$Rdiff["upper"],
                                           maxlower=outWald2$Rdiff["lower"],                        
                                           whichAlgo=contr$algo)
        #--
        outEmpLikeSR <- EmpLikSRCIKM(data=d,
                                     tstar=t,
                                     mytol=contr$tol,
                                     alpha=1-level,
                                     maxupper=outWald1S$RR["upper"],
                                     minlower=outWald1S$RR["lower"],
                                     minupper=outWald2S$RR["upper"],
                                     maxlower=outWald2S$RR["lower"],                        
                                     whichAlgo=contr$algo)
        # }}}
        # {{{ Fix possible NA in EL CI (if search interval is not large enough)
        # {{{ RR
        PbRRl <- is.na(outEmpLikeRR[1])
        PbRRu <- is.na(outEmpLikeRR[2])
        PbRRboth <- all(c(PbRRl,PbRRu))
        PbRR <- any(c(PbRRl,PbRRu))
        if(PbRR){
            while(PbRR & k < 10){           
                k <- k+1
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current is:  ",paste(format(outEmpLikeRR,digits=5),collapse=" - ")))
                }
                #-- recompute limits to search for EL CI, based on Wald CI  -- #
                outWald1k <- ciRisk(R1=risk1,
                                    R2=risk0,
                                    seR1=se.risk1,
                                    seR2=se.risk0,
                                    RR.H0=RR.H0,
                                    Diff.H0=Diff.H0,
                                    alpha=(1-level)/(k**1.3))    
                outWald2k <- ciRisk(R1=risk1,
                                    R2=risk0,
                                    seR1=se.risk1,
                                    seR2=se.risk0,
                                    RR.H0=RR.H0,
                                    Diff.H0=Diff.H0,
                                    alpha=(1-level)*k)                    
                #-
                if(contr$Trace){
                    if(PbRRboth) print(paste0("Now search limits are, for lower= ",paste(format(c(outWald1k$RR["lower"],outWald2k$RR["lower"]),digits=5),collapse=" - "),", for  upper= ",paste(format(c(outWald2k$RR["upper"],outWald1k$RR["upper"]),digits=5),collapse=" - ")))
                    if(PbRRl) print(paste0("Now search limits are lower= ",paste(format(c(outWald1k$RR["lower"],outWald2k$RR["lower"]),digits=5),collapse=" - ")))
                    if(PbRRu) print(paste0("Now search limits are upper= ",paste(format(c(outWald2k$RR["upper"],outWald1k$RR["upper"]),digits=5),collapse=" - ")))
                }
                #-- recompute EL CI  -- #
                outEmpLikeRR <- EmpLikRRCIKM(data=d,
                                             tstar=t,
                                             mytol=contr$tol,
                                             alpha=1-level,
                                             maxupper=ifelse(PbRRu,outWald1k$RR["upper"],outEmpLikeRR[2]+100*contr$tol),
                                             minlower=ifelse(PbRRl,outWald1k$RR["lower"],outEmpLikeRR[1]-100*contr$tol),
                                             minupper=ifelse(PbRRu,outWald2k$RR["upper"],outEmpLikeRR[2]-100*contr$tol),
                                             maxlower=ifelse(PbRRl,outWald2k$RR["lower"],outEmpLikeRR[1]+100*contr$tol),
                                             whichAlgo=contr$algo)
                if(contr$Trace){
                    print(paste0("Updated is:  ",paste(format(outEmpLikeRR,digits=5),collapse=" - ")))
                }
                PbRRl <- is.na(outEmpLikeRR[1])
                PbRRu <- is.na(outEmpLikeRR[2])
                PbRRboth <- all(c(PbRRl,PbRRu))
                PbRR <- any(c(PbRRl,PbRRu))
            }
        }
        # }}}  
        # {{{ Rdiff
        k <- contr$k
        PbRdiffl <- is.na(outEmpLikeRdiff[1])
        PbRdiffu <- is.na(outEmpLikeRdiff[2])
        PbRdiffboth <- all(c(PbRdiffl,PbRdiffu))
        PbRdiff <- any(c(PbRdiffl,PbRdiffu))    
        if(PbRdiff){
            while(PbRdiff & k < 10){           
                k <- k+1
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current is:  ",paste(format(outEmpLikeRdiff,digits=5),collapse=" - ")))
                }
                #-- recompute limits to search for EL CI, based on Wald CI  -- #
                outWald1k <- ciRisk(R1=risk1,
                                    R2=risk0,
                                    seR1=se.risk1,
                                    seR2=se.risk0,
                                    RR.H0=RR.H0,
                                    Diff.H0=Diff.H0,
                                    alpha=(1-level)/(k**1.3))    
                outWald2k <- ciRisk(R1=risk1,
                                    R2=risk0,
                                    seR1=se.risk1,
                                    seR2=se.risk0,
                                    RR.H0=RR.H0,
                                    Diff.H0=Diff.H0,
                                    alpha=(1-level)*k)                    
                #-
                if(contr$Trace){
                    if(PbRdiffboth) print(paste0("Now search limits are, for lower= ",paste(format(c(outWald1k$Rdiff["lower"],outWald2k$Rdiff["lower"]),digits=5),collapse=" - "),", for  upper= ",paste(format(c(outWald2k$Rdiff["upper"],outWald1k$Rdiff["upper"]),digits=5),collapse=" - ")))
                    if(PbRdiffl) print(paste0("Now search limits are lower= ",paste(format(c(outWald1k$Rdiff["lower"],outWald2k$Rdiff["lower"]),digits=5),collapse=" - ")))
                    if(PbRdiffu) print(paste0("Now search limits are upper= ",paste(format(c(outWald2k$Rdiff["upper"],outWald1k$Rdiff["upper"]),digits=5),collapse=" - ")))
                }
                #-- recompute EL CI  -- #
                outEmpLikeRdiff <- EmpLikRdiffCIKM(data=d,
                                                   tstar=t,
                                                   mytol=contr$tol,
                                                   alpha=1-level,
                                                   maxupper=ifelse(PbRdiffu,outWald1k$Rdiff["upper"],outEmpLikeRdiff[2]+100*contr$tol),
                                                   minlower=ifelse(PbRdiffl,outWald1k$Rdiff["lower"],outEmpLikeRdiff[1]-100*contr$tol),
                                                   minupper=ifelse(PbRdiffu,outWald2k$Rdiff["upper"],outEmpLikeRdiff[2]-100*contr$tol),
                                                   maxlower=ifelse(PbRdiffl,outWald2k$Rdiff["lower"],outEmpLikeRdiff[1]+100*contr$tol),
                                                   whichAlgo=contr$algo)
                if(contr$Trace){
                    print(paste0("Updated is:  ",paste(format(outEmpLikeRdiff,digits=5),collapse=" - ")))
                }
                ## print(outEmpLikeRdiff)
                PbRdiffl <- is.na(outEmpLikeRdiff[1])
                PbRdiffu <- is.na(outEmpLikeRdiff[2])
                PbRdiffboth <- all(c(PbRdiffl,PbRdiffu))
                PbRdiff <- any(c(PbRdiffl,PbRdiffu))
            }
        }
        # }}}

        # {{{ SR
        k <- contr$k
        PbSRl <- is.na(outEmpLikeSR[1])
        PbSRu <- is.na(outEmpLikeSR[2])
        PbSRboth <- all(c(PbSRl,PbSRu))
        PbSR <- any(c(PbSRl,PbSRu))
        #--
        if(PbSR){
            while(PbSR & k < 10){       
                k <- k+1
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current CI for SR is:  ",paste(format(outEmpLikeSR,digits=5),collapse=" - ")))
                }
                #-- recompute limits to search for EL CI, based on Wald CI  -- #
                outWald1Sk <- ciRisk(R1=surv1,
                                     R2=surv0,
                                     seR1=se.risk1,
                                     seR2=se.risk0,
                                     RR.H0=RR.H0,
                                     Diff.H0=Diff.H0,
                                     alpha=(1-level)/(k**1.3))    
                outWald2Sk <- ciRisk(R1=surv1,
                                     R2=surv0,
                                     seR1=se.risk1,
                                     seR2=se.risk0,
                                     RR.H0=RR.H0,
                                     Diff.H0=Diff.H0,
                                     alpha=(1-level)*k)                    
                #-
                if(contr$Trace){
                    if(PbSRboth){
                        print(paste0("Now search limits are, for lower= ",paste(format(c(outWald1Sk$RR["lower"],outWald2Sk$RR["lower"]),digits=5),collapse=" - "),", for  upper= ",paste(format(c(outWald2Sk$RR["upper"],outWald1Sk$RR["upper"]),digits=5),collapse=" - ")))
                    }else{
                        if(PbSRl) print(paste0("Now search limits are lower= ",paste(format(c(outWald1Sk$RR["lower"],outWald2Sk$RR["lower"]),digits=5),collapse=" - ")))
                        if(PbSRu) print(paste0("Now search limits are upper= ",paste(format(c(outWald2Sk$RR["upper"],outWald1Sk$RR["upper"]),digits=5),collapse=" - ")))
                    }
                }
                #-- recompute EL CI  -- #
                outEmpLikeSR <- EmpLikSRCIKM(data=d,
                                             tstar=t,
                                             mytol=contr$tol,
                                             alpha=1-level,
                                             maxupper=ifelse(PbSRu,outWald1Sk$RR["upper"],outEmpLikeSR[2]+100*contr$tol),
                                             minlower=ifelse(PbSRl,outWald1Sk$RR["lower"],outEmpLikeSR[1]-100*contr$tol),
                                             minupper=ifelse(PbSRu,outWald2Sk$RR["upper"],outEmpLikeSR[2]-100*contr$tol),
                                             maxlower=ifelse(PbSRl,outWald2Sk$RR["lower"],outEmpLikeSR[1]+100*contr$tol),
                                             whichAlgo=contr$algo)
                if(contr$Trace){
                    print(paste0("Updated is:  ",paste(format(outEmpLikeSR,digits=5),collapse=" - ")))
                }
                #--
                PbSRl <- is.na(outEmpLikeSR[1])
                PbSRu <- is.na(outEmpLikeSR[2])
                PbSRboth <- all(c(PbSRl,PbSRu))
                PbSR <- any(c(PbSRl,PbSRu))
            }
        }
        # }}}
   
        # }}}
        # }}}

        # {{{ EL CI p-values
        p.Diff.EL <- LikeRatioRRRdiffKM(data=d,Rdiff=Diff.H0,tstar=t,mytol=contr$tol,whichAlgo=contr$algo)$pval
        p.RR.EL <- LikeRatioRRRdiffKM(data=d,RR=RR.H0,tstar=t,mytol=contr$tol,whichAlgo=contr$algo)$pval
        p.SR.EL <- LikeRatioRRRdiffKM(data=d,SR=SR.H0,tstar=t,mytol=contr$tol,whichAlgo=contr$algo)$pval
        # }}}

    }else{
        p.Diff.EL <- p.RR.EL <- p.SR.EL <- NA
        outEmpLikeSR <- outEmpLikeRR <- outEmpLikeRdiff <- c(NA,NA)
    }
    # {{{ output
    table.SR <- rbind(Wald=res.Wald.SR$RR[c(1,3:5)],
                      EL=c(res.Wald.SR$RR["est."],outEmpLikeSR,p.SR.EL))
    table.RR <- rbind(Wald=res.Wald$RR[c(1,3:5)],
                      EL=c(res.Wald$RR["est."],outEmpLikeRR,p.RR.EL))
    table.Diff <- rbind(Wald=res.Wald$Rdiff[c(1,3:5)],
                        EL=c(res.Wald$Rdiff["est."],outEmpLikeRdiff,p.Diff.EL))
    colnames(table.SR)[4] <- colnames(table.Diff)[4] <- colnames(table.RR)[4] <- "p-value"
    #--
    x <- list(table.RR=table.RR,
              table.Diff=table.Diff,
              table.SR=table.SR,
              input=list(t=t,
                         RR.H0=RR.H0,
                         Diff.H0=Diff.H0,
                         SR.H0=SR.H0,
                         level=level,                         
                         contr=contr)
              )
    class(x) <- "TwoSampleKaplanMeier"
    # }}}
    x 
}




#----------------------------------------------------------------------
### TwoSampleKaplanMeier.R ends here
