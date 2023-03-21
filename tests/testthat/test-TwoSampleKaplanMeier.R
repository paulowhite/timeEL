### test-TwoSampleKaplanMeier.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:41) 
## Version: 
## Last-Updated: Aug 12 2022 (13:13) 
##           By: Paul Blanche
##     Update #: 6
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
    #--
    Res2SKM95 <- TwoSampleKaplanMeier(time=Freireich$time,
                                      status=Freireich$status,
                                      group=Freireich$group,
                                      t=10)
    Res2SKM90 <- TwoSampleKaplanMeier(time=Freireich$time,
                                      status=Freireich$status,
                                      group=Freireich$group,
                                      level=0.9,
                                      t=10)
    # See Table 4 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(Res2SKM95$table.SR["EL",c("est.","lower","upper")], c(1.976,1.139,3.973),tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(Res2SKM90$table.SR["EL",c("est.","lower","upper")], c(1.976,1.241,3.498),tolerance=5e-4,ignore_attr = TRUE)
})



test_that("Freireich data, time=20 weeks: Thomas and Grunkemeier (JASA, 1975)",{
    #--
    Res2SKM95 <- TwoSampleKaplanMeier(time=Freireich$time,
                                      status=Freireich$status,
                                      group=Freireich$group,
                                      t=20)
    Res2SKM90 <- TwoSampleKaplanMeier(time=Freireich$time,
                                      status=Freireich$status,
                                      group=Freireich$group,
                                      level=0.9,
                                      t=20)
    # See Table 4 in Thomas and Grunkemeier (JASA, 1975)
    expect_equal(Res2SKM95$table.SR["EL",c("est.","lower","upper")], c(6.591,2.145,39.02),tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(Res2SKM90$table.SR["EL",c("est.","lower","upper")], c(6.591,2.507,27.41),tolerance=5e-4,ignore_attr = TRUE)
})


# }}}



#----------------------------------------------------------------------
### test-TwoSampleKaplanMeier.R ends here
