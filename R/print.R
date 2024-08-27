# {{{ print.KaplanMeier
##' @title Print function for object of class 'KaplanMeier'
##' 
##' @param x an object of class 'KaplanMeier'
##' @param digits number of digits to print the results
##' @param type either "surv" or "risk" (the default), depending on whether we want to print the results in terms of a survival probability or a risk (i.e., one minus the survival probability).
##' @param method either "EL", "Wald" or "both", depending on whether we want to print the results obtained when using empirical likelihood inference (EL), Wald-type inference (Wald) or both. Default is 'NULL', which means that 'method' inherits the value of the corresponding control parameter used when creating the object 'x'. 
##' @param ... Not used
##' 
##' @author Paul Blanche
##'
##' @examples
##' # This example reproduces some results presented in Table 1 of Thomas and Grunkemeier (1975)
##' ResKM.1.95 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
##'                           status=Freireich$status[Freireich$group==1],
##'                           t=10, level=0.95, contr=list(tol=1e-4))
##' print(ResKM.1.95, digits=3, type="surv", method="EL")   # EL results for survival
##' print(ResKM.1.95, digits=3, type="risk", method="Wald") # Wald results for risk
##'
##' @return no return value, called for printing only.
##' 
##' @export
print.KaplanMeier <- function(x,digits=4,
                              type="risk",
                              method=NULL,...){
    #--
    if(is.null(method)){
        method <- x$input$contr$method
    }
    #--
    if( method=="EL" & x$input$contr$method=="Wald"){
        stop("No results are available to print for method 'EL'.")
    }
    #--
    WhichWord <- "risk"
    #--
    if(type=="surv"){
        x$table[,c("risk","lower","upper")] <- 1-x$table[,c("risk","upper","lower")]
        WhichWord <- "survival probability"
        if(!is.null(x$input$risk.H0)){
            x$input$risk.H0 <- 1-x$input$risk.H0
        }
    }
    #---
    cat("\n")
    cat("Non-parametric inference on the ",WhichWord," at time t=",specdec(x$input$t,digits),". \n")
    cat("\n")
    cat("--------")
    if(method=="Wald" | method=="both"){
        cat("\n")
        cat("Wald-type inference (via cloglog transformation):\n")
        cat("\n")
        cat("Estimate=",specdec(x$table["Wald","risk"],digits)," and ",
            100*x$input$level, "% confidence interval= [",specdec(x$table["Wald","lower"],digits),";",specdec(x$table["Wald","upper"],digits),"] \n")    
        if(!is.null(x$input$risk.H0)){
            if(x$table["Wald","p-value"]>10^(-digits)){
                thep <- specdec(x$table["Wald","p-value"],digits)
            }else{
                thep <- paste0("<", 10^(-digits))
            }
            cat("\n")
            cat("p-value=",thep," for 'H0: ",WhichWord,"=",specdec(x$input$risk.H0,digits),"'.")
            cat("\n")
        }
    }
    if(method=="both"){
        cat("\n")
        cat("--------")
    }
    if(method=="EL" | method=="both"){
        cat("\n")
        cat("Empirical likelihood inference:\n")
        cat("\n")
        cat("Estimate=",specdec(x$table["EL","risk"],digits)," and ",
            100*x$input$level, "% confidence interval= [",specdec(x$table["EL","lower"],digits),";",specdec(x$table["EL","upper"],digits),"] \n")
        if(!is.null(x$input$risk.H0)){
            if(x$table["EL","p-value"]>10^(-digits)){
                thep <- specdec(x$table["EL","p-value"],digits)
            }else{
                thep <- paste0("<", 10^(-digits))
            }
            cat("\n")
            cat("p-value=",thep," for 'H0: ",WhichWord,"=",specdec(x$input$risk.H0,digits),"'.")
            cat("\n")
        }
        cat("\n")
    }
    return(invisible(NULL))
}
# }}}


# {{{ print.TwoSampleKaplanMeier
##' @title Print function for object of class 'TwoSampleKaplanMeier'
##' 
##' @param x an object of class 'TwoSampleKaplanMeier'
##' @param digits number of digits to print the results
##' @param what either "SR", "RR", "Diff" or "all" (default), depending on whether we want to print the results for the survival ratio (SR), the risk ratio (RR), the risk difference (Diff) or all of them.
##' @param method either "EL", "Wald" or "both", depending on whether we want to print the results obtained when using empirical likelihood inference (EL), Wald-type inference (Wald) or both. Default is 'NULL', which means that 'method' inherits the value of the corresponding control parameter used when creating the object 'x'. 
##' @param ... Not used
##' 
##' @author Paul Blanche
##'
##' @examples
##' # This example reproduces some results presented in Table 4 of Thomas and Grunkemeier (1975)
##' Res2SKM95 <- TwoSampleKaplanMeier(time=Freireich$time,
##'                                   status=Freireich$status,
##'                                   group=Freireich$group,
##'                                   t=10)
##' print(Res2SKM95, digits=3, what="SR", method="EL")
##'
##' @return no return value, called for printing only.
##' 
##' @export
print.TwoSampleKaplanMeier <- function(x,
                                       digits=4,
                                       what="all", # "all", "SR", "RR" or "Diff"
                                       method=NULL,...){# "both", "EL" or "Wald"
    #--
    if(is.null(method)){
        method <- x$input$contr$method
    }
    #--
    if( method=="EL" & x$input$contr$method=="Wald"){
        stop("No results are available to print for method 'EL'.")
    }
    #--
    # {{{ RR, Diff
    if(what %in% c("all","RR","Diff")){
        y <- x
        class(y) <- "TwoSampleAalenJohansen"
        what2 <- what
        if(what=="all"){
            what2 <- "both"
        }
        print(y,digits,what2,method,absRisk=FALSE)
    }
    # }}}
    if(what=="all"){
        cat("--------")
        cat("\n")
    }
    # {{{ SR
    if(what %in% c("all","SR")){
        cat("\n")
        if(method=="Wald" | method=="both"){
            cat("Wald-type inference for the Survival Ratio (via log transformation):\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.SR["Wald","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.SR["Wald","lower"],digits),";",specdec(x$table.SR["Wald","upper"],digits),"] \n")    
            if(!is.null(x$input$SR.H0)){
                if(x$table.SR["Wald","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.SR["Wald","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: SR=",specdec(x$input$SR.H0,digits),"'.")
                cat("\n")
            }
        }
        if(method=="both"){
            cat("\n")
            cat("--------")
            cat("\n")
        }
        if(method=="EL" | method=="both"){
            cat("Empirical likelihood inference for the Survival Ratio:\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.SR["EL","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.SR["EL","lower"],digits),";",specdec(x$table.SR["EL","upper"],digits),"] \n")    
            if(!is.null(x$input$SR.H0)){
                if(x$table.SR["EL","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.SR["EL","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: SR=",specdec(x$input$SR.H0,digits),"'.")
                cat("\n")
            }    
            cat("\n")
        }
    }
    # }}}
    return(invisible(NULL))
}
# }}}
    

# {{{ print.AalenJohansen
##' @title Print function for object of class 'AalenJohansen'
##' 
##' @param x an object of class 'AalenJohansen'
##' @param digits number of digits to print the results
##' @param method either "EL", "Wald" or "both", depending on whether we want to print the results obtained when using empirical likelihood inference (EL), Wald-type inference (Wald) or both. Default is 'NULL', which means that 'method' inherits the value of the corresponding control parameter used when creating the object 'x'. 
##' @param ... Not used
##' 
##' @author Paul Blanche
##'
##' @examples
##' x <- AalenJohansen(time=melanoma5$time, cause=melanoma5$status, t=4, level=0.95)
##' print(x, digits=3, method="EL")
##'
##' @return no return value, called for printing only.
##' 
##' @export
print.AalenJohansen <- function(x,
                                digits=4,
                                method=NULL,...){
    #--
    if(is.null(method)){
        method <- x$input$contr$method
    }
    #--
    if( method=="EL" & x$input$contr$method=="Wald"){
        stop("No results are available to print for method 'EL'.")
    }
    #--
    cat("\n")
    cat("Non-parametric inference on the (absolute) risk at time t=",specdec(x$input$t,digits),". \n")
    cat("\n")
    cat("--------")
    if(method=="Wald" | method=="both"){
        cat("\n")
        cat("Wald-type inference (via cloglog transformation):\n")
        cat("\n")
        cat("Estimate=",specdec(x$table["Wald","risk"],digits)," and ",
            100*x$input$level, "% confidence interval= [",specdec(x$table["Wald","lower"],digits),";",specdec(x$table["Wald","upper"],digits),"] \n")    
        if(!is.null(x$input$risk.H0)){
            if(x$table["Wald","p-value"]>10^(-digits)){
                thep <- specdec(x$table["Wald","p-value"],digits)
            }else{
                thep <- paste0("<", 10^(-digits))
            }
            cat("\n")
            cat("p-value=",thep," for 'H0: risk=",specdec(x$input$risk.H0,digits),"'.")
            cat("\n")
        }
    }
    if(method=="both"){
        cat("\n")
        cat("--------")
        cat("\n")
    }
    if(method=="EL" | method=="both"){
        cat("Empirical likelihood inference:\n")
        cat("\n")
        cat("Estimate=",specdec(x$table["EL","risk"],digits)," and ",
            100*x$input$level, "% confidence interval= [",specdec(x$table["EL","lower"],digits),";",specdec(x$table["EL","upper"],digits),"] \n")    
        if(!is.null(x$input$risk.H0)){
            if(x$table["EL","p-value"]>10^(-digits)){
                thep <- specdec(x$table["EL","p-value"],digits)
            }else{
                thep <- paste0("<", 10^(-digits))
            }
            cat("\n")
            cat("p-value=",thep," for 'H0: risk=",specdec(x$input$risk.H0,digits),"'.")
            cat("\n")
        }
        cat("\n")
    }
    return(invisible(NULL))
}
# }}}


# {{{ print.TwoSampleAalenJohansen
##' @title Print function for object of class 'TwoSampleAalenJohansen'
##' 
##' @param x an object of class 'TwoSampleAalenJohansen'
##' @param digits number of digits to print the results
##' @param what either "RR", "Diff" or "both" (default), depending on whether we want to print the results for the risk ratio (RR), the risk difference (Diff) or both.
##' @param method either "EL", "Wald" or "both", depending on whether we want to print the results obtained when using empirical likelihood inference (EL), Wald-type inference (Wald) or both. Default is 'NULL', which means that 'method' inherits the value of the corresponding control parameter used when creating the object 'x'. 
##' @param absRisk Default is TRUE and this should not be changed.
##' @param ... Not used
##' 
##' @author Paul Blanche
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
##' print(ResSimA100, digits=3, what="Diff")
##' print(ResSimA100, digits=3, what="RR")
##' 
##' @examplesIf FALSE
##' ## Same example data, but now analyzed with and empirical likelihood approach. It
##' ## takes approx 20 seconds to run.
##' 
##' ResSimA100 <- TwoSampleAalenJohansen(time=SimA100$time,
##'                                      cause=SimA100$status,
##'                                      group=SimA100$group,
##'                                      t=1)
##' print(ResSimA100, digits=3, what="Diff",  method="EL")
##'
##' @return no return value, called for printing only.
##' 
##' @export
print.TwoSampleAalenJohansen <- function(x,
                                         digits=4,
                                         what="both", # "both", "RR" or "Diff"
                                         method=NULL, # "both", "EL" or "Wald"
                                         absRisk=TRUE,...){ # just to better deal with this function being used also in the non-competing risk case
    #--
    if(is.null(method)){
        method <- x$input$contr$method
    }
    #--
    if( method=="EL" & x$input$contr$method=="Wald"){
        stop("No results are available to print for method 'EL'.")
    }
    #--
    cat("\n")
    if(absRisk){
        cat("Non-parametric inference to compare two (absolute) risks at time t=",specdec(x$input$t,digits),". \n")
    }else{
        cat("Non-parametric inference to compare two risks (or survival probabilities) at time t=",specdec(x$input$t,digits),". \n")
    }
    cat("\n")
    cat("--------")
    # {{{ Risk Ratio
    if(what=="both" | what=="RR"){
        cat("\n")
        if(method=="Wald" | method=="both"){
            cat("Wald-type inference for the Risk Ratio (via log transformation):\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.RR["Wald","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.RR["Wald","lower"],digits),";",specdec(x$table.RR["Wald","upper"],digits),"] \n")    
            if(!is.null(x$input$RR.H0)){
                if(x$table.RR["Wald","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.RR["Wald","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: RR=",specdec(x$input$RR.H0,digits),"'.")
                cat("\n")
            }
        }
        if(method=="both"){
            cat("\n")
            cat("--------")
            cat("\n")
        }
        if(method=="EL" | method=="both"){
            cat("Empirical likelihood inference for the Risk Ratio:\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.RR["EL","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.RR["EL","lower"],digits),";",specdec(x$table.RR["EL","upper"],digits),"] \n")    
            if(!is.null(x$input$RR.H0)){
                if(x$table.RR["EL","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.RR["EL","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: RR=",specdec(x$input$RR.H0,digits),"'.")
                cat("\n")
            }    
            cat("\n")
        }
    }
    # }}}
    if(what=="both"){
        cat("--------")
        cat("\n")
    }
    # {{{ Risk difference
    if(what=="both" | what=="Diff"){
        cat("\n")
        if(method=="Wald" | method=="both"){
            cat("Wald-type inference for the Risk Difference:\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.Diff["Wald","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.Diff["Wald","lower"],digits),";",specdec(x$table.Diff["Wald","upper"],digits),"] \n")    
            if(!is.null(x$input$Diff.H0)){
                if(x$table.Diff["Wald","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.Diff["Wald","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: Diff=",specdec(x$input$Diff.H0,digits),"'.")
                cat("\n")
            }
        }
        if(method=="both"){
            cat("\n")
            cat("--------")
            cat("\n")
        }
        if(method=="EL" | method=="both"){
            cat("Empirical likelihood inference for the Risk Difference:\n")
            cat("\n")
            cat("Estimate=",specdec(x$table.Diff["EL","est."],digits)," and ",
                100*x$input$level, "% confidence interval= [",specdec(x$table.Diff["EL","lower"],digits),";",specdec(x$table.Diff["EL","upper"],digits),"] \n")    
            if(!is.null(x$input$Diff.H0)){
                if(x$table.Diff["EL","p-value"]>10^(-digits)){
                    thep <- specdec(x$table.Diff["EL","p-value"],digits)
                }else{
                    thep <- paste0("<", 10^(-digits))
                }
                cat("\n")
                cat("p-value=",thep," for 'H0: Diff=",specdec(x$input$Diff.H0,digits),"'.")
                cat("\n")
            }    
            cat("\n")
        }
    }
    # }}}
    return(invisible(NULL))
}
# }}}
