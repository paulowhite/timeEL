### MyLogLikeKM.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 31 2021 (13:13) 
## Version: 
## Last-Updated: Dec 21 2021 (11:45) 
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


# This works only with no ties...
## MyLogLikeKM <- function(a,ni){
    ## sum(log(a)) + sum((ni-1)*log(1 - a))
## }

# This works also with ties
MyLogLikeKM <- function(a,ni,di){
    sum(di*log(a)) + sum((ni-di)*log(1 - a))
}

#----------------------------------------------------------------------
### MyLogLikeKM.R ends here
