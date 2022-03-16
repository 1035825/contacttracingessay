t <- 0.3 #fixed tlat

fixed_iter_delay = function(d,B){ #p_c* as a function of delta, R0
  n = floor((1-d)/(t-d)) #n as defined in essay
  control = function(p){ #1 - new reproduction number R as a function of p
    1 - (1-p)*B*((p*B)^n - 1)/(p*B-1) #looking for this to equal 0 (R=1)
  }
  intercept = uniroot(control,c(0,1.5)) 
  return(unname(unlist(intercept[1]))) #p_c* such that R=1
}

X4 <- seq(0.001,0.299,by=0.001) 
Y4 <- rep.int(0,299) 
for(i in 1:299){
  Y4[i] <- fixed_iter_delay(X4[i],2) #p_c* for R0=2, 0.001<delta<0.299
}
DF4 <- data.frame(X4,Y4)
X3 <- seq(0.001,0.299,by=0.001) 
Y3 <- rep.int(0,299) 
for(i in 1:299){
  Y3[i] <- fixed_iter_delay(X3[i],1.3) #p_c* for R0=1.3, 0.001<delta<0.299
}
DF3 <- data.frame(X3,Y3)
X10 <- seq(0.001,0.299,by=0.001) 
Y10 <- rep.int(0,299) 
for(i in 1:299){
  Y10[i] <- fixed_iter_delay(X10[i],1.1) #p_c* for R0=1.1, 0.001<delta<0.299
}
DF10 <- data.frame(X10,Y10)

pdf("iteration2.pdf", width=6, height=6)
plot(DF4,xlim=c(0,0.33),ylim=c(0,1),type="l",lwd=2,xlab="Tracing delay",ylab="Critical tracing efficiency",col="red")
par(new=TRUE)
plot(DF3,xlim=c(0,0.33),ylim=c(0,1),type="l",lwd=2,xlab="Tracing delay",ylab="Critical tracing efficiency",col="blue")
par(new=TRUE)
plot(DF10,xlim=c(0,0.33),ylim=c(0,1),type="l",lwd=2,xlab="Tracing delay",ylab="Critical tracing efficiency",col="black")
abline(v=0.3,lty=2)
legend(0.01,0.15, legend=c("R0 = 1.1", "R0 = 1.3", "R0 = 2"),col=c("black", "blue", "red"), lty=1:1:1, cex=0.8)
dev.off()








