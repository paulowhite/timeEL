### test-AalenJohansen.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:42) 
## Version: 
## Last-Updated: Aug 12 2022 (15:10) 
##           By: Paul Blanche
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ check AalenJohansen function results equivalent to KaplanMeier function results when no competing risks
test_that("check AalenJohansen function results equivalent to KaplanMeier function results when no competing risks",{

    AJres <- AalenJohansen(time=melanoma5$time,
                           cause=as.numeric(melanoma5$status!=0)*(3-melanoma5$status),
                           t=0.1,
                           level=0.95,contr=list(method="EL",tol=1e-7))


    AJresKM <- KaplanMeier(time=melanoma5$time,
                           status=as.numeric(melanoma5$status!=0),
                           t=0.1,
                           level=0.95,contr=list(method="both",tol=1e-7,Trace=TRUE))
    ## max(abs(AJresKM$table-AJres$table),na.rm=TRUE)   
    expect_equal(AJres$table,
                 AJresKM$table,
                 tolerance=1e-6) 
})
# }}}


# {{{ Test Aalen-Johansen routines: reproduce results of Blanche (LIDA, 2020)
test_that("Melanoma data: Blanche (LIDA, 2020)",{
    alltime <- c(1.5,2,3,3.5,4,5,6,8,10,12)
    nt <- length(alltime)
    EL.upper <- EL.lower <- Wald.upper <- Wald.lower <- rep(NA,nt)
    for(thet in alltime){
        AJres <- AalenJohansen(time=melanoma5$time,
                               cause=melanoma5$status,
                               t=thet,
                               level=0.95)
        Wald.lower[which(thet==alltime)] <- AJres$table["Wald","lower"]
        Wald.upper[which(thet==alltime)] <- AJres$table["Wald","upper"]
        EL.lower[which(thet==alltime)] <- AJres$table["EL","lower"]
        EL.upper[which(thet==alltime)] <- AJres$table["EL","upper"]        
    }
    # Check Table 2 in Blanche (LIDA, 2020)    
    expect_equal(Wald.lower*100, c(1.2,1.6,7.2,9.1,9.6,13.2,16.9,21.5,24.1,24.1),tolerance=5e-3)
    expect_equal(Wald.upper*100, c(6.8,7.6,16.8,19.4,20.1,25.0,30.5,37.4,42.1,42.1),tolerance=5e-3)
    expect_equal(EL.lower*100, c(1.0,1.4,6.9,8.8,9.3,12.9,16.5,21.1,23.8,23.8),tolerance=5e-2) # A bug (n instead of n-1) has been fixed since he publication, hence a minor difference, and the tolerance set to tolerance=5e-2.
    expect_equal(EL.upper*100, c(6.1,6.9,16.2,18.9,19.5,24.4,29.9,36.8,41.6,41.6),tolerance=5e-2) # A bug (n instead of n-1) has been fixed since he publication, hence a minor difference, and the tolerance set to tolerance=5e-2.
})
# }}}

#----------------------------------------------------------------------
### test-AalenJohansen.R ends here
