finaldelayalter=function(R){ #same code as in section 8.3, except setting delta=tlat=0 everywhere
  #so exp(d) <- 1, etc.
  long <- function(p){
    P_0 = function(x){
      prod = rep.int(0,x+1)
      for(i in 1:(x+1)){
        prod[i] <- i+p*R
      }
      prodfinal <- 1
      for(i in 1:(x+1)){
        prodfinal <- prodfinal*prod[i]
      }
      ((p*R)^x)*(x+1)/prodfinal
    }
    P_1 = function(x){
      prod = rep.int(0,x+1)
      for(i in 1:(x+1)){
        prod[i] <- 1+i+p*R
      }
      prodfinal <- 1
      for(i in 1:(x+1)){
        prodfinal <- prodfinal*prod[i]
      }
      ((p*R)^x)*(x+2)/prodfinal
    }
    vect <- rep.int(0,100)
    for(x in 1:100){
      vect[x] <- (x)*P_0(x)
    }
    k_10 <- sum(vect)
    k_00 <- ((1-p)/p)*k_10
    VECT <- rep.int(0,100)
    for(x in 1:100){
      VECT[x] <- (x)*P_1(x)
    }
    k_11 <- sum(VECT)
    k_01 <- ((1-p)/p)*k_11
    
    value <- function(r){
      (k_00 - r)*(k_11 - r) - k_01*k_10
    }
    RR <- uniroot(value,c(0.1,10))
    sol <- unname(unlist(RR[1]))
    return(sol)
  }
  
  Rvect <- rep.int(0,1000)
  for(i in 1:1000){
    Rvect[i] <- long(i/1000)
  }
  p_c <- (which(abs(Rvect - 1) == min(abs(Rvect - 1))))/1000
  return(p_c)
}

Rvalues <- rep.int(0,100)
for(i in 1:100){
  Rvalues[i] <- finaldelayalter(1.1+i/100) #p_c* for single-step tracing, 1.11<R0<2.10
}
Rvalues[Rvalues == 1] <- NA #removing any entries where p_c*>1

iterative <- function(x){
  1-1/x #calculating p_c* for iterative tracing, as defined in essay
}
x <- rep.int(0,100)
for(i in 1:100){
  x[i] <- iterative(1.1+i/100) #p_c* for iterative tracing, 1.11<R0<2.10
}

pdf("B3.pdf",height=6,width=6)
plot(seq(1.11,2.10,by=0.01),Rvalues[1:100],type="l",ylim=c(0,1),col="red",xlab="R0",ylab="Critical tracing efficiency")
par(new=TRUE)
plot(seq(1.11,2.10,by=0.01),x[1:100],type="l",ylim=c(0,1),col="black",xlab="R0",ylab="Critical tracing efficiency")
legend(1.15,0.95,legend=c("Single-step", "Iterative"),col=c("red","black"), lty=1:1, cex=0.8)
dev.off()
