### test-AJwithSE.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:39) 
## Version: 
## Last-Updated: Aug 12 2022 (14:58) 
##           By: Paul Blanche
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ Check internal CIF and se computation against that of prodlim
test_that("Check internal CIF and se computation against that of prodlim",{
    thet <- 15
    AJmine <- AJwithSE(tstar=thet,data=melanoma5,CompAllse=TRUE)
    AJprodlim <- prodlim::prodlim(prodlim::Hist(time,status)~1,melanoma5)
    #--
    MyRisk <- AJmine$CIF[-1]
    MyseRisk <- AJmine$se[-1]
    ProdlimRisk <- AJprodlim$cuminc$`1`
    ProdlimRisk <- ProdlimRisk[AJprodlim$time<=thet]
    ProdlimseRisk <- AJprodlim$se.cuminc$`1`
    ProdlimseRisk <- ProdlimseRisk[AJprodlim$time<=thet]
    #-- check CIF
    expect_equal(MyRisk,ProdlimRisk,tolerance=5e-10,ignore_attr = TRUE)
    # check se
    expect_equal(ProdlimseRisk,MyseRisk,ignore_attr = TRUE,tolerance=5e-10) 
})
# }}}

#----------------------------------------------------------------------
### test-AJwithSE.R ends here
