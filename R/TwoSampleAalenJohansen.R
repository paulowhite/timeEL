### TwoSampleAalenJohansen.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 15 2021 (09:30) 
## Version: 
## Last-Updated: Sep 13 2023 (11:39) 
##           By: Paul Blanche
##     Update #: 287
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

##' @description
##' Computes an (absolute) risk difference or ratio with right-censored competing risks data, together with a confidence interval and a
##' p-value (to test for a difference between the two risks). Pointwise estimates are computed via the Aalen-Johansen estimator. Computation of confidence intervals
##' and p-values are based on either Empirical Likelihood (EL) inference or Wald-type inference. Both are non-parametric approaches, which are asymptotically equivalent.
##' For the Wald-type approach, the asymptotic normal approximation is used on the log scale for the risk ratio. No transformation is used for the risk difference.
##' See Blanche & Eriksson (2023) for details.
##'
##' 
##' 
##' @title Risk difference and ratio using the Aalen-Johansen method
##' 
##' @param time vector of times (possibly censored)
##' @param cause vector of event types/causes. It should be coded 1 for main events, 2 for competing events and 0 for censored.
##' @param group vector of binary group indicator. The reference group should be coded 0, the other 1. 
##' @param t the time point of interest (e.g. 1 to compute a 1-year risk ratio)
##' @param RR.H0 the risk ratio under the null hypothesis, to compute a p-value. Default is 1.
##' @param Diff.H0 the risk difference under the null hypothesis, to compute a p-value. Default is 0.
##' @param level confidence level for the confidence intervals. Default is 0.95. 
##' @param contr list of control parameters. tol=tolerance for numerical computation, default is 1e-5. method="EL", "Wald" or "both"  indicates wether 95% CI and the p-value should be computed based on Empirical Likelihood (EL) inference, Wald-type inference or both. algo=2 (default) or 1, depending on which computational method should be used to maximize the empirical likelihood (method 1 or 2, as described in Blanche & Eriksson (2023))
##'
##' @return an object of class 'TwoSampleAalenJohansen'
##' 
##' @author Paul Blanche
##'
##' @references
##' Blanche & Eriksson (2023). Empirical likelihood comparison of absolute risks.
##'
##' @examples
##' ## A simple example for Wald-type inference, using simulated data.
##' ## It illustrates the possible inconsistency of Wald-type inference, in
##' ## terms of statistical significance, when inference is based on the risk
##' ## ratio and on the risk difference. This inconsistency cannot exist
##' ## using an empirical likelihood approach.
##' 
##' ResSimA100 <- TwoSampleAalenJohansen(time=SimA100$time,
##'                                      cause=SimA100$status,
##'                                      group=SimA100$group,
##'                                      t=1,
##'                                      contr=list(method="Wald"))
##' ResSimA100
##' 
##' @examplesIf FALSE
##' ## Same example data, but now analyzed with and empirical likelihood approach. It
##' ## takes approx 20 seconds to run.
##' 
##' ResSimA100 <- TwoSampleAalenJohansen(time=SimA100$time,
##'                                      cause=SimA100$status,
##'                                      group=SimA100$group,
##'                                      t=1)
##' ResSimA100
##'
##'
##' @export
TwoSampleAalenJohansen <- function(time,
                                   cause,
                                   group,
                                   t,
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
    if(!("numeric" %in% class(time))){
        stop("time should be numeric.")
    }
    if(!all(cause %in% c(0,1,2))){
        stop("cause should contain values which are either 0, 1 or 2.")
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
    if(!(Diff.H0 > -1 & Diff.H0 < 1)){
        stop("Diff.H0 should be within ]-1,1[.")
    }
    if( !(level > 0 & level < 1)   ){
        stop("level should be within ]0,1[.")
    }
    # }}}
    
    # {{{ data 
    d <- data.frame(time=time,status=cause,group=group)
    d <- d[order(d$group,d$time),]
    
    # {{{ To handle ties in times
    if(contr$method  %in% c("both","EL")){
        #---
        d0 <- d[d$group==0,]
        d1 <- d[d$group==1,]
        # Do for each group
        for(j in 1:2){
            dddd  <- list(d0,d1)[[j]]
            # check for ties in times 
            if(any(duplicated(dddd$time))){
                WhereAreTies <- which(duplicated(dddd$time))
                WhereMax <- max(which(dddd$time<=t))
                # check for ties in times before time t
                if(any(WhereAreTies <= WhereMax)){
                    # if any before t, problem. We stop.
                    stop("The function, used with method='EL' or method='both', cannot handle ties in times, when they occur before t. You can use the function with method='Wald', though.")
                }else{
                    # if all after t, no problem.
                    # As the values of the times after t do not matter, as long as they are after t, we randomly update the data to break the ties.
                    dddd$time[WhereAreTies] <- max(dddd$time)+stats::runif(length(WhereAreTies))
                    dddd <- dddd[order(dddd$time),]                
                }
            }
            # update data
            if(j==1) d0 <- dddd
            if(j==2) d1 <- dddd            
        }
        # merge updated data
        d <- rbind.data.frame(d0,d1)
        d <- d[order(d$group,d$time),]
    }
    # }}}     
    # check at least one event in each group
    if(min(d$time[d$status==1 & d$group==0])>t | min(d$time[d$status==1 & d$group==1])>t){
        stop("There is no main event observed in the data before time t in at least one of the two groups.")
    }    
    # }}}
  
    # {{{ Wald inference
    res.AJwithSE.0 <- AJwithSE(tstar=t,data=d[d$group==0,],CompAllse=FALSE)
    res.AJwithSE.1 <- AJwithSE(tstar=t,data=d[d$group==1,],CompAllse=FALSE) 
    risk1 <- res.AJwithSE.1$risk.t
    risk0 <- res.AJwithSE.0$risk.t
    se.risk1 <- res.AJwithSE.1$se.t
    se.risk0 <- res.AJwithSE.0$se.t
    #---
    res.Wald <- ciRisk(R1=risk1,
                       R2=risk0,
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
    
        outWald2 <- ciRisk(R1=risk1,
                           R2=risk0,
                           seR1=se.risk1,
                           seR2=se.risk0,
                           RR.H0=RR.H0,
                           Diff.H0=Diff.H0,
                           alpha=(1-level)*k)
        # }}}
        # {{{ EL CI inference
        # {{{ Computation with first search interval
        outEmpLikeRR <- EmpLikRRCI(data=d,
                                   tstar=t,
                                   mytol=contr$tol,
                                   alpha=1-level,
                                   maxupper=outWald1$RR["upper"],
                                   minlower=outWald1$RR["lower"],
                                   minupper=outWald2$RR["upper"],
                                   maxlower=outWald2$RR["lower"],                        
                                   whichAlgo=contr$algo)    
        outEmpLikeRdiff <- EmpLikRdiffCI(data=d,
                                         tstar=t,
                                         mytol=contr$tol,
                                         alpha=1-level,
                                         maxupper=outWald1$Rdiff["upper"],
                                         minlower=outWald1$Rdiff["lower"],
                                         minupper=outWald2$Rdiff["upper"],
                                         maxlower=outWald2$Rdiff["lower"],                        
                                         whichAlgo=contr$algo)
        # }}}
        # {{{ Fix possible NA in EL CI (if search interval is not large enough)
        # {{{ RR
        PbRRl <- is.na(outEmpLikeRR[1])
        PbRRu <- is.na(outEmpLikeRR[2])
        PbRRboth <- all(c(PbRRl,PbRRu))
        PbRR <- any(c(PbRRl,PbRRu))
        ## browser()
        if(PbRR){
            while(PbRR & k < 11){           
                k <- k+1
                if(k==10){k <- 100}
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current CI for RR is:  ",paste(format(outEmpLikeRR,digits=5),collapse=" - ")))
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
                                    alpha=min((1-level)*k,0.9))                    
                #-
                if(contr$Trace){
                    if(PbRRboth) print(paste0("Now search limits are (for CI of RR) for lower= ",paste(format(c(outWald1k$RR["lower"],outWald2k$RR["lower"]),digits=5),collapse=" - "),", for  upper= ",paste(format(c(outWald2k$RR["upper"],outWald1k$RR["upper"]),digits=5),collapse=" - ")))
                    if(PbRRl) print(paste0("Now search limits are (for CI of RR) lower= ",paste(format(c(outWald1k$RR["lower"],outWald2k$RR["lower"]),digits=5),collapse=" - ")))
                    if(PbRRu) print(paste0("Now search limits are (for CI of RR) upper= ",paste(format(c(outWald2k$RR["upper"],outWald1k$RR["upper"]),digits=5),collapse=" - ")))
                }
                #-- recompute EL CI  -- #
                outEmpLikeRR <- EmpLikRRCI(data=d,
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
                #--
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
            while(PbRdiff & k < 11){           
                k <- k+1
                if(k==10){k <- 100}
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current CI for Diff is:  ",paste(format(outEmpLikeRdiff,digits=5),collapse=" - ")))
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
                                    alpha=min((1-level)*k,0.9))                    
                #-
                if(contr$Trace){
                    if(PbRdiffboth) print(paste0("Now search limits are (for CI of Diff) for lower= ",paste(format(c(outWald1k$Rdiff["lower"],outWald2k$Rdiff["lower"]),digits=5),collapse=" - "),", for  upper= ",paste(format(c(outWald2k$Rdiff["upper"],outWald1k$Rdiff["upper"]),digits=5),collapse=" - ")))
                    if(PbRdiffl) print(paste0("Now search limits are (for CI of Diff) lower= ",paste(format(c(outWald1k$Rdiff["lower"],outWald2k$Rdiff["lower"]),digits=5),collapse=" - ")))
                    if(PbRdiffu) print(paste0("Now search limits are (for CI of Diff) upper= ",paste(format(c(outWald2k$Rdiff["upper"],outWald1k$Rdiff["upper"]),digits=5),collapse=" - ")))
                }
                #-- recompute EL CI  -- #
                outEmpLikeRdiff <- EmpLikRdiffCI(data=d,
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
        # }}}
        # }}}
        # {{{ EL CI p-values
        p.Diff.EL <- LikeRatioRRRdiff(data=d,Rdiff=Diff.H0,tstar=t,mytol=contr$tol,whichAlgo=contr$algo)$pval
        p.RR.EL <- LikeRatioRRRdiff(data=d,RR=RR.H0,tstar=t,mytol=contr$tol,whichAlgo=contr$algo)$pval
        # }}}    
    }else{
        p.Diff.EL <- p.RR.EL <- NA
        outEmpLikeRR <- outEmpLikeRdiff <- c(NA,NA)
    }
    # {{{ output
    table.RR <- rbind(Wald=res.Wald$RR[c(1,3:5)],
                      EL=c(res.Wald$RR["est."],outEmpLikeRR,p.RR.EL))
    table.Diff <- rbind(Wald=res.Wald$Rdiff[c(1,3:5)],
                        EL=c(res.Wald$Rdiff["est."],outEmpLikeRdiff,p.Diff.EL))
    colnames(table.Diff)[4] <- colnames(table.RR)[4] <- "p-value"
    #--
    x <- list(table.RR=table.RR,
              table.Diff=table.Diff,
              input=list(t=t,
                         RR.H0=RR.H0,
                         Diff.H0=Diff.H0,
                         level=level,                         
                         contr=contr)
              )
    class(x) <- "TwoSampleAalenJohansen"
    # }}}
    x 
}


#----------------------------------------------------------------------
### TwoSampleAalenJohansen.R ends here
