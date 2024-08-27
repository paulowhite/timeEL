plot.constrainedAJ <- function(x, plot.lambda0=TRUE){
    k <- length(x$time)
    xx <- c(0,rep(x$time,each=2),x$tstar)
    # CIF 1
    yAJcons <- c(rep(c(0,x$CIF1c[-k]),each=2),x$CIF1c[k],x$CIF1star)
    yAJ0 <- x$outCIFlambda0$CIF1
    yAJ <-  rep(c(0,yAJ0),each=2)
    # CIF 2
    yAJcons2 <- rep(c(0,x$CIF2c),each=2)
    yAJ02 <- x$outCIFlambda0$CIF2
    yAJ2 <-  rep(c(0,yAJ02),each=2)
    # CIF 1    
    graphics::plot.default(xx,yAJcons,type="l",col="red",ylim=c(0,1),xlim=c(0,1.1*x$tstar),lwd=2,
                           xlab="time", ylab="Cumulative Incidence function")
    if(plot.lambda0){
        graphics::lines(xx,yAJ,col="red",lwd=2,lty=2)
    }
    # CIF 2
    graphics::lines(xx,yAJcons2,col="blue",lwd=2,lty=1)
    if(plot.lambda0){
        graphics::lines(xx,yAJ2,col="blue",lwd=2,lty=2)
    }
    # legend
    graphics::abline(h=x$CIF1star,lty=2)
    graphics::abline(v=x$tstar,lty=2)
    graphics::legend("topleft",
                     col=rep(c("red","blue"),each=2),
                     legend=c("AJ constrained cause=1","AJ cause=1",
                              "AJ constrained cause=2","AJ cause=2"),
                     lty=rep(1:2,times=2),
                     lwd=2)
}

plot.constrainedKM <- function(x){
    # fix if one value for Sc
    k <- length(x$time)
    xx <-c(0,rep(x$time,each=2),x$tstar)
    yKMcons <- c(rep(c(1,x$Sc[-k]),each=2),x$Sc[k],x$pstar)
    yKM0 <- x$x$KM[x$x$time <= x$tstar]
    yKM <-  rep(c(1,yKM0),each=2)
    cbind(xx, yKMcons, yKM)
    graphics::plot.default(xx,yKMcons,type="l",col="red",ylim=c(0,1),xlim=c(0,1.1*x$tstar),lwd=2)
    graphics::lines(xx,yKM,col="blue",lwd=2)
    graphics::abline(h=x$pstar,lty=2)
    graphics::abline(v=x$tstar,lty=2)
    graphics::legend("bottomleft",col=c("red","blue"),legend=c("KM","KM constrained"),lty=1,lwd=2)
}


plot.AJlambda <- function(x){
    k <- length(x$time)
    xx <- c(0,rep(x$time,each=2),x$tstar)
    # CIF 1
    yAJcons <- rep(c(0,x$CIF1c),each=2)
    # CIF 2
    yAJcons2 <- rep(c(0,x$CIF2c),each=2)
    # CIF 1    
    graphics::plot.default(xx,yAJcons,type="l",col="red",ylim=c(0,1),xlim=c(0,1.1*x$tstar),lwd=2,
                           xlab="time", ylab="Cumulative Incidence function")
    # CIF 2
    graphics::lines(xx,yAJcons2,col="blue",lwd=2,lty=1)
    # legend
    graphics::abline(v=x$tstar,lty=2)
    graphics::abline(h=x$CIFstar,lty=2)
    graphics::legend("topleft",
                     col=c("red","blue"),
                     legend=c("AJ constrained cause=1",
                              "AJ constrained cause=2"),
                     lwd=2, bty="n",
                     title=paste("CIFstar=",round(x$CIFstar,3),", lambda=",round(x$lambda,5)))
}
