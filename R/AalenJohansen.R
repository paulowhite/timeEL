### AalenJohansen.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 12 2021 (16:11) 
## Version: 
## Last-Updated: Aug 10 2023 (13:56) 
##           By: Paul Blanche
##     Update #: 163
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

##' @description
##' Computes the Aalen-Johansen estimator to estimate an (absolute) risk
##' with right-censored competing risks data, together with a confidence interval and (possibly) a
##' p-value (for a one-sample hypothesis test). Computation of confidence intervals
##' and p-values is based on either Empirical Likelihood (EL) inference
##' or Wald-type inference. Both are non-parametric approaches, which are asymptotically equivalent.
##' See Blanche (2020) for details. For the Wald-type approach, the asymptotic normal approximation is used on the cloglog scale. See e.g. equation 4.21 in Beyersmann et al (2011).
##' 
##' @title Risk estimate using the Aalen-Johansen method
##' @param time vector of times (possibly censored)
##' @param cause vector of event types/causes. It should be coded 1 for main events, 2 for competing events and 0 for censored.
##' @param t the time point of interest (e.g., 1 to compute a 1-year risk)
##' @param risk.H0 risk under the null hypothesis, if one would like to compute the correspondng p-value. Default is NULL, for which no p-value is computed.
##' @param level confidence level for the confidence intervals. Default is 0.95. 
##' @param contr list of control parameters. tol=tolerance for numerical computation, default is 1e-5. method="EL", "Wald" or "both"  indicates wether 95% CI and p-value should be computed based on Empirical Likelihood (EL) inference , Wald-type inference or both.
##' @return object of class 'AalenJohansen'
##' @author Paul Blanche
##'
##' @examples
##' x <- AalenJohansen(time=melanoma5$time, cause=melanoma5$status, t=4, level=0.95)
##' x
##'
##' @references
##' Blanche, P. (2020). Confidence intervals for the cumulative incidence function via constrained NPMLE. Lifetime Data Analysis, 26(1), 45-64.
##'
##' Beyersmann, Allignol,  & Schumacher (2011). Competing risks and multistate models with R. Springer Science & Business Media.
##' 
##' @export
AalenJohansen <- function(time,
                          cause,
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
    if(!all(cause %in% c(0,1,2))){
        stop("cause should contain values which are either 0, 1 or 2.")
    }
    if(t <=0 ){
        stop("The time point of interest, t, should be positive.")
    }
    if(!is.null(risk.H0)){
        if( !(risk.H0 > 0 & risk.H0 < 1)){
            stop("risk.H0 should be within ]0,1[.")
        }
    }
    # }}}
    
    # data 
    d <- data.frame(time=time,status=cause)
    d <- d[order(d$time),]

    # {{{ To handle ties in times
    if(contr$method  %in% c("both","EL")){
        # check for ties in times 
        if(any(duplicated(d$time))){
            WhereAreTies <- which(duplicated(d$time))
            WhereMax <- max(which(d$time<=t))
            # check for ties in times before time t
            if(any(WhereAreTies <= WhereMax)){
                # if any before t, problem. We stop.
                stop("The function, used with method='EL' or method='both', cannot handle ties in times, when they occur before t. You can use the function with method='Wald', though.")
            }else{
                # if all after t, no problem.
                # As the values of the times after t do not matter, as long as they are after t, we randomly update the data to break the ties.
                d$time[WhereAreTies] <- max(d$time)+stats::runif(length(WhereAreTies))
                d <- d[order(d$time),]                
            }
        }
    }
    ## browser()
    # }}}
    
    # check at least one event
    if(min(d$time[d$status==1])>t){
        stop("There is no main event observed in the data before time t.")
    }
    # Wald CI  (cloglog scale)
    AJfit <- AJwithSE(tstar=t,data=d,CompAllse=FALSE)
    risk <- AJfit$risk.t
    se.risk <- AJfit$se.t   
    CI.Wald <- c(lower=lowerRisk(est=risk,se=se.risk,alpha=1-level),upper=upperRisk(est=risk,se=se.risk,alpha=1-level))
    se.clolog.scale <- -se.risk/((1-risk)*log(1-risk))
    est.clolog.scale <- log(-log(1-risk))

    # EL CI
    if(contr$method=="both" | contr$method=="EL"){
        k <- contr$k       
        CI.EL <- LRCIAJ(t,
                        d,
                        alpha=1-level,
                        maxupper=upperRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3)), # maximum for grid search of upper limit
                        minlower=lowerRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3)), # minimum for grid search of lower limit
                        minupper=upperRisk(est=risk,se=se.risk,alpha=(1-level)*k), # minimum for grid search of upper limit
                        maxlower=lowerRisk(est=risk,se=se.risk,alpha=(1-level)*k),
                        mytol=contr$tol)
        # {{{ Fix possible NA in EL CI (if search interval is not large enough)
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
                if(Pbu & !Pbl){
                    Newmaxupper <- upperRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newminupper <- upperRisk(est=risk,se=se.risk,alpha=(1-level)*k)
                    Newminlower <- CI.EL[1]-100*contr$tol
                    Newmaxlower <- CI.EL[1]+100*contr$tol
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits for  upper= ",paste(format(c(Newminupper,Newmaxupper),digits=5),collapse=" - ")))
                    }
                }
                if(Pbl & !Pbu){
                    Newmaxupper <- CI.EL[2]+100*contr$tol
                    Newminupper <- CI.EL[2]-100*contr$tol
                    Newminlower <- lowerRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newmaxlower <- lowerRisk(est=risk,se=se.risk,alpha=(1-level)*k)
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits for lower= ",paste(format(c(Newminlower,Newmaxlower),digits=5),collapse=" - ")))
                    }
                }
                if(Pbboth){
                    Newmaxupper <- upperRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newminupper <- upperRisk(est=risk,se=se.risk,alpha=(1-level)*k)
                    Newminlower <- lowerRisk(est=risk,se=se.risk,alpha=(1-level)/(k**1.3))
                    Newmaxlower <- lowerRisk(est=risk,se=se.risk,alpha=(1-level)*k)
                    #--
                    if(contr$Trace){
                        print(paste0("Now search limits are, for lower= ",paste(format(c(Newminlower,Newmaxlower),digits=5),collapse=" - "),
                                     ", for  upper= ",paste(format(c(Newminupper,Newmaxupper),digits=5),collapse=" - ")))
                    }
                }
                #-- recompute limits to search for EL CI, based on Wald CI  -- #
                CI.EL <- LRCIAJ(t,
                                d,
                                alpha=1-level,
                                maxupper=Newmaxupper,
                                minlower=Newminlower,
                                minupper=Newminupper,
                                maxlower=Newmaxlower,
                                mytol=contr$tol)
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
            p.EL <- LR.test.AJ(t, risk.H0,  d)$p
        }else{
            p.EL <- NA
        }
    }else{
        p.EL <- p.Wald <-NA
    }
    # output
    ## browser()
    res.Wald <- c(risk=risk,CI.Wald,"p-value"=p.Wald)
    res.EL <- c(risk=risk,CI.EL,"p-value"=p.EL)
    table <- rbind(Wald=res.Wald,EL=res.EL)
    x <- list(table=table,
              input=list(risk.H0=risk.H0,
                         t=t,
                         level=level,
                         contr=contr)
              )
    class(x) <- "AalenJohansen"
    x 
}

#----------------------------------------------------------------------
### AalenJohansen.R ends here
