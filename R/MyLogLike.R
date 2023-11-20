### MyLogLike.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Dec  9 2020 (12:53) 
## Version: 
## Last-Updated: Aug 14 2023 (09:50) 
##           By: Paul Blanche
##     Update #: 20
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
    # {{{ fix, just to handle the case in which the survival goes to 0
    # and to avoid computing 0*log(0), which returns NaN,
    # although it should return 0 here.
    ToSumSurvTerm <- (ssn:(ssn-nbt+1))*log(1 - a1 - a2)
    if(nbt>ssn){ ToSumSurvTerm <- ToSumSurvTerm[-nbt] }
    # }}}    
    # loglikelihood corresponding to AJ1
    ll1 <- sum(log(a1[i1])) + sum(log(a2[i2])) + sum(ToSumSurvTerm)
    ll1
}


#----------------------------------------------------------------------
### MyLogLike.R ends here
