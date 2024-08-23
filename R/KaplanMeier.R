### KaplanMeier.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 30 2021 (14:53) 
## Version: 
## Last-Updated: Aug  1 2024 (15:28) 
##           By: Paul Blanche
##     Update #: 155
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @description
##' Computes the Kaplan-Meier estimator to estimate a risk or, equivalently, a survival probability,
##' with right-censored data, together with a confidence interval and (possibly) a
##' p-value (for a one-sample hypothesis test). Computation of confidence intervals
##' and p-values is based on either Empirical Likelihood (EL) inference
##' or Wald-type inference. Both are non-parametric approaches, which are asymptotically equivalent.
##' See Thomas & Grunkemeier (1975) for details about the Empirical Likelihood method. For the Wald-type approach,
##' the asymptotic normal approximation is used on the cloglog scale. See e.g. equation 4.16 in Beyersmann et al (2011).
##' 
##' @title Risk and survival probability estimates using the Kaplan-Meier method
##' 
##' @param time vector of times (possibly censored)
##' @param status vector of usual survival status indicators (0 for censored observations, 1 for events)
##' @param t the time point of interest (e.g. 1 to compute a 1-year risk or survival probability)
##' @param risk.H0 risk under the null hypothesis, if one would like to compute the correspondng p-value. Default is NULL, for which no p-value is computed.
##' @param level confidence level for the confidence intervals. Default is 0.95. 
##' @param contr list of control parameters. tol=tolerance for numerical computation, default is 1e-5. method="EL", "Wald" or "both"  indicates wether 95% CI and p-value should be computed based on Empirical Likelihood (EL) inference , Wald-type inference or both.
##' @return object of class 'KaplanMeier'
##' @author Paul Blanche
##'
##' @references
##' Thomas & Grunkemeier (1975). Confidence interval estimation of survival probabilities for censored data. Journal of the American Statistical Association, 70(352), 865-871.
##' 
##' Beyersmann, Allignol,  & Schumacher (2011). Competing risks and multistate models with R. Springer Science & Business Media.
##' @examples
##' # This example reproduces some results presented in Table 1 of Thomas and Grunkemeier (1975)
##' ResKM.1.95 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
##'                           status=Freireich$status[Freireich$group==1],
##'                           t=10, level=0.95, contr=list(tol=1e-4))
##' ResKM.1.95
##'
##' @export
KaplanMeier <- function(time,
                        status,
                        t,
                        risk.H0=NULL,
                        level=0.95,
                        contr=list(tol=1e-5,k=3,Trace=FALSE,method="both")){


    # {{{ About input parameters
    # to handle the fact that we might want to give only some of the arguments in contr
    if(!missing(contr)){
        if(is.null(contr$tol))  contr$tol <- 1e-5       
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
    if(!all(status %in% c(0,1))){
        stop("status should contain values which are either 0 or 1.")
    }
    if(t <=0 ){
        stop("The time point of interest, t, should be positive.")
    }
    if(!is.null(risk.H0)){
        if( !(risk.H0 > 0 & risk.H0 < 1)   ){
            stop("risk.H0 should be within ]0,1[.")
        }
    }
    # }}}
   
    # data
    d <- data.frame(time=time,status=status)
    # check at least one event
    if(min(d$time[d$status==1])>t){
        stop("There is no event observed in the data before time t.")
    }
    # KM estimates
    KMfit <- KMwithSE(tstar=t,data=d)
    risk <- 1-KMfit$res.tstar$Surv
    # Greewood s.e.
    se.risk <- KMfit$res.tstar$se  
    # Wald CI  (cloglog scale)
    CI.Wald <- c(lower=lowerRisk(est=risk,se=se.risk,alpha=1-level),upper=upperRisk(est=risk,se=se.risk,alpha=1-level))
    se.clolog.scale <- -se.risk/((1-risk)*log(1-risk))
    est.clolog.scale <- log(-log(1-risk))
    # EL CI
    ## browser()
    if(contr$method=="both" | contr$method=="EL"){
        k <- contr$k       
        CI.EL <- rev(1-LRCIKM(tstar=t,
                              data=d,
                              mytol=contr$tol,
                              alpha=1-level,
                              maxupper=upperRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3)), # maximum for grid search of upper limit
                              minlower=lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3)), # minimum for grid search of lower limit
                              minupper=upperRisk(est=1-risk,se=se.risk,alpha=(1-level)*k), # minimum for grid search of upper limit
                              maxlower=lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)*k)
                              ))
        # {{{ Fix possible NA in EL CI (if search interval is not large enough)
        # # RK: upper and lower in 'maxupper', 'minupper' etc refer to CI of  survival, not CI of risk....
        Pbl <- is.na(CI.EL[1])
        Pbu <- is.na(CI.EL[2])
        Pbboth <- all(c(Pbl,Pbu))
        Pb <- any(c(Pbl,Pbu))
        ## browser()
        if(Pb){
            while(Pb & k < 10){           
                k <- k+1
                if(contr$Trace){
                    print(paste0("Wile loop, k=",k))
                    print(paste0("Current CI is:  ",paste(format(CI.EL,digits=5),collapse=" - ")))
                }
                if(!Pbu & Pbl){ 
                    Newmaxupper <- upperRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newminupper <- upperRisk(est=1-risk,se=se.risk,alpha=(1-level)*k)
                    Newminlower <- (1-CI.EL[2])-100*contr$tol
                    Newmaxlower <- (1-CI.EL[2])+100*contr$tol
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits for  upper= ",paste(format(c(Newminupper,Newmaxupper),digits=5),collapse=" - ")))
                    }
                }
                if(!Pbl & Pbu){
                    Newmaxupper <- (1-CI.EL[1])+100*contr$tol
                    Newminupper <- (1-CI.EL[1])-100*contr$tol
                    Newminlower <- lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newmaxlower <- lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)*k)
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits for lower= ",paste(format(c(Newminlower,Newmaxlower),digits=5),collapse=" - ")))
                    }
                }
                if(Pbboth){
                    Newmaxupper <- upperRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newminupper <- upperRisk(est=1-risk,se=se.risk,alpha=(1-level)*k)
                    Newminlower <- lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newmaxlower <- lowerRisk(est=1-risk,se=se.risk,alpha=(1-level)*k)
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits are, for lower= ",paste(format(c(Newminlower,Newmaxlower),digits=5),collapse=" - "),
                                     ", for  upper= ",paste(format(c(Newminupper,Newmaxupper),digits=5),collapse=" - ")))
                    }
                }
                #-- recompute limits to search for EL CI, based on Wald CI  -- #
                CI.EL <- rev(1-LRCIKM(tstar=t,
                                      data=d,
                                      mytol=contr$tol,
                                      alpha=1-level,
                                      maxupper=Newmaxupper,
                                      minlower=Newminlower,
                                      minupper=Newminupper,
                                      maxlower=Newmaxlower
                                      ))
                # Is there still a problem?
                Pbl <- is.na(CI.EL[1])
                Pbu <- is.na(CI.EL[2])
                Pbboth <- all(c(Pbl,Pbu))
                Pb <- any(c(Pbl,Pbu))
            }
        }

        # }}}
    }else{
        CI.EL <- c(NA,NA)
    }
    # p-value if needed
    if(!is.null(risk.H0)){
        p.Wald <- 2*stats::pnorm(abs(est.clolog.scale-log(-log(1-risk.H0)))/se.clolog.scale,lower.tail=FALSE)
        if(contr$method=="both" | contr$method=="EL"){
            p.EL <- LR.test.KM(tstar=t, pstar=1-risk.H0,  data=d)$p
        }else{
            p.EL <- NA
        }
    }else{
        p.EL <- p.Wald <-NA
    }
    # output
    res.Wald <- c(risk=risk,CI.Wald,"p-value"=p.Wald)
    res.EL <- c(risk=risk,CI.EL,"p-value"=p.EL)
    table <- rbind(Wald=res.Wald,EL=res.EL)
    x <- list(table=table,
              input=list(risk.H0=risk.H0,
                         t=t,
                         level=level,
                         contr=contr)
              )
    class(x) <- "KaplanMeier"
    x 
}

#----------------------------------------------------------------------
### KaplanMeier.R ends here
