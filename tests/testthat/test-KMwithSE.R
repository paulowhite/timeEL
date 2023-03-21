### test-KMwithSE.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:40) 
## Version: 
## Last-Updated: Aug 12 2022 (12:08) 
##           By: Paul Blanche
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ Check the internal Greenwood variance computation of KM against that of prodlim
test_that("Check the internal Greenwood variance computation of KM against that of prodlim",{

    #-- test with Freireich data-------
    ResKM.1.gw <- KMwithSE(tstar=10,data=Freireich)
    ResKM.1.prodlim <- prodlim::prodlim(prodlim::Hist(time,status)~1,Freireich)
    expect_equal(ResKM.1.gw$res.full$KM,
                 ResKM.1.prodlim$surv,
                 tolerance=5e-10)
    expect_equal(sqrt(ResKM.1.gw$res.full$gw),
                 ResKM.1.prodlim$se.surv,
                 tolerance=5e-10) 
    #-- test with simulated data -------
    dsim <- prodlim::SimSurv(3000)
    ResKM.dsim.gw <- KMwithSE(tstar=10,data=dsim)
    ResKM.dsim.prodlim <- prodlim::prodlim(prodlim::Hist(time,status)~1,dsim)
    #--
    expect_equal(ResKM.dsim.gw$res.full$KM,
                 ResKM.dsim.prodlim$surv,
                 tolerance=5e-10)
    # note, we need to remove the results from the last time (prodlim creates some NaN)
    expect_equal(sqrt(ResKM.dsim.gw$res.full$gw)[-length(ResKM.dsim.prodlim$se.surv)],
                 ResKM.dsim.prodlim$se.surv[-length(ResKM.dsim.prodlim$se.surv)],
                 tolerance=5e-10) 
})

# }}}

#----------------------------------------------------------------------
### test-KMwithSE.R ends here
