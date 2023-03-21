### MyLogLike.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec  9 2020 (12:53) 
## Version: 
## Last-Updated: Aug 11 2022 (14:36) 
##           By: Paul Blanche
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

MyLogLike <- function(a1,
                      a2,
                      ssn
                      ){
    i1 <- which(a1 !=0)
    i2 <- which(a2 !=0)
    # number of observerved times before the timepoint of ineterest
    nbt <- length(a1)
    ssn <- ssn-1
    # loglikelihood corresponding to AJ1
    ll1 <- sum(log(a1[i1])) + sum(log(a2[i2])) + sum((ssn:(ssn-nbt+1))*log(1 - a1 - a2))
    ll1
}


#----------------------------------------------------------------------
### MyLogLike.R ends here
