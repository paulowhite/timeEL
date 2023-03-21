### test-KaplanMeier.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:40) 
## Version: 
## Last-Updated: Sep  6 2022 (10:22) 
##           By: Paul Blanche
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:




# {{{ Test Kaplan-Meier routines: reproduce results of Thomas and Grunkemeier (JASA, 1975)

test_that("Freireich data, time=10 weeks: Thomas and Grunkemeier (JASA, 1975)",{  
    # Active group
    ResKM.1.95 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
                              status=Freireich$status[Freireich$group==1],
                              t=10,
                              level=0.95,
                              contr=list(tol=1e-4))
    ResKM.1.90 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
                              status=Freireich$status[Freireich$group==1],
                              t=10,
                              level=0.90,
                              contr=list(tol=1e-4))
    # See Table 1 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(ResKM.1.95$table["EL",c("risk","lower","upper")], 1-c(0.753,0.904,0.540),tolerance=6e-4,ignore_attr = TRUE) 
    expect_equal(ResKM.1.90$table["EL",c("risk","lower","upper")], 1-c(0.753,0.885,0.576),tolerance=7e-4,ignore_attr = TRUE)
    # control group
    ResKM.0.95 <- KaplanMeier(time=Freireich$time[Freireich$group==0],
                              status=Freireich$status[Freireich$group==0],
                              t=10,
                              level=0.95,
                              contr=list(tol=1e-4))
    ResKM.0.90 <- KaplanMeier(time=Freireich$time[Freireich$group==0],
                              status=Freireich$status[Freireich$group==0],
                              t=10,
                              level=0.90,
                              contr=list(tol=1e-4))
    # See Table 1 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(ResKM.0.95$table["EL",c("risk","lower","upper")], 1-c(0.381,0.593,0.196),tolerance=5e-4,ignore_attr = TRUE)  
    expect_equal(ResKM.0.90$table["EL",c("risk","lower","upper")], 1-c(0.381,0.560,0.222),tolerance=5e-4,ignore_attr = TRUE) 
})


test_that("Freireich data, time=20 weeks: Thomas and Grunkemeier (JASA, 1975)",{
    # Active group
    ResKM.1.95 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
                              status=Freireich$status[Freireich$group==1],
                              t=20,
                              level=0.95,
                              contr=list(tol=1e-4))
    ResKM.1.90 <- KaplanMeier(time=Freireich$time[Freireich$group==1],
                              status=Freireich$status[Freireich$group==1],
                              t=20,
                              level=0.90,
                              contr=list(tol=1e-4))
    # See Table 1 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(ResKM.1.95$table["EL",c("risk","lower","upper")], 1-c(0.628,0.821,0.395),tolerance=5e-3,ignore_attr = TRUE)  
    expect_equal(ResKM.1.90$table["EL",c("risk","lower","upper")], 1-c(0.628,0.795,0.432),tolerance=5e-3,ignore_attr = TRUE) 
    # control group
    ResKM.0.95 <- KaplanMeier(time=Freireich$time[Freireich$group==0],
                              status=Freireich$status[Freireich$group==0],
                              t=20,
                              level=0.95,
                              contr=list(tol=1e-4))
    ResKM.0.90 <- KaplanMeier(time=Freireich$time[Freireich$group==0],
                              status=Freireich$status[Freireich$group==0],
                              t=20,
                              level=0.90,
                              contr=list(tol=1e-4))
    # See Table 1 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(ResKM.0.95$table["EL",c("risk","lower","upper")], 1-c(0.095,0.266,0.016),tolerance=5e-4,ignore_attr = TRUE) 
    expect_equal(ResKM.0.90$table["EL",c("risk","lower","upper")], 1-c(0.095,0.233,0.023),tolerance=5e-4,ignore_attr = TRUE) 
})
# }}}

# {{{ compare with package km.ci
test_that("Compare to km.ci results",{
    dsim <- prodlim::SimSurv(50)
    # Active group: my function
    ResKM.1.95 <- KaplanMeier(time=dsim$time,
                              status=dsim$status,
                              t=4,
                              level=0.95,
                              contr=list(tol=1e-4))
    # Active group: my function
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data=dsim)
    fit2 <- km.ci::km.ci(fit,method="grunkemeier")
    reskmci <- 1-c(summary(fit2,time=4)$surv,summary(fit2,time=4)$upper,summary(fit2,time=4)$lower)
    # See Table 1 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(ResKM.1.95$table["EL",c("risk","lower","upper")],
                 reskmci,
                 tolerance=5e-3,ignore_attr = TRUE)      
})
# }}}



#----------------------------------------------------------------------
### test-KaplanMeier.R ends here
