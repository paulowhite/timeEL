### test-TwoSampleAalenJohansen.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 11 2022 (16:41) 
## Version: 
## Last-Updated: Aug 10 2023 (13:10) 
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

# {{{ Example with simulated data (Scenario A, n=100) and Wald-type inference inconsistency

test_that("Reproduce example with simulated data (Scenario A, n=100) and Wald-type inference inconsistency",{

    ResSimA100 <- TwoSampleAalenJohansen(time=SimA100$time,
                                         cause=SimA100$status,
                                         group=SimA100$group,
                                         t=1)

    # Check Wald-type inference results    
    expect_equal(c(round(ResSimA100$table.RR["Wald",1:3],2),round(ResSimA100$table.RR["Wald",4],3)),c(2.52,0.93,6.82,0.070) ,tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(c(round(ResSimA100$table.Diff["Wald",1:3],2),round(ResSimA100$table.Diff["Wald",4],3)),c(0.24,0.03,0.44,0.026) ,tolerance=5e-4,ignore_attr = TRUE)
    # Check EL inference results
    expect_equal(c(round(ResSimA100$table.RR["EL",1:3],2),round(ResSimA100$table.RR["EL",4],3)),c(2.52,1.04,8.17,0.039) ,tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(c(round(ResSimA100$table.Diff["EL",1:3],2),round(ResSimA100$table.Diff["EL",4],3)),c(0.24,0.01,0.43,0.039) ,tolerance=5e-4,ignore_attr = TRUE)
})
# }}}



# {{{ Reproduce example with BMTtcell data (see Figure in Blanche & Eriksson (2023))

test_that("Reproduce example with BMTtcell data and Wald-type inference inconsistency",{

    Restcell1 <- TwoSampleAalenJohansen(time=BMTtcell$time,
                                        cause=BMTtcell$status,
                                        group=BMTtcell$group,
                                        t=15)

    # Check Wald-type inference results    
    expect_equal(c(round(Restcell1$table.RR["Wald",1:3],2),round(Restcell1$table.RR["Wald",4],3)),c(0.630, 0.380, 1.020, 0.062) ,tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(c(round(Restcell1$table.Diff["Wald",1:3],2),round(Restcell1$table.Diff["Wald",4],3)),c(-0.150, -0.270, -0.020, 0.024) ,tolerance=5e-4,ignore_attr = TRUE)
    # Check EL inference results
    expect_equal(c(round(Restcell1$table.RR["EL",1:3],2),round(Restcell1$table.RR["EL",4],3)),c(0.630, 0.360, 0.970, 0.036) ,tolerance=5e-4,ignore_attr = TRUE)
    expect_equal(c(round(Restcell1$table.Diff["EL",1:3],2),round(Restcell1$table.Diff["EL",4],3)),c(-0.150, -0.260, -0.010,  0.036) ,tolerance=5e-4,ignore_attr = TRUE)
})
# }}}





#----------------------------------------------------------------------
### test-TwoSampleAalenJohansen.R ends here
