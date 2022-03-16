R <- 1.5

final=function(t){
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
  vect <- rep.int(0,100)
  for(x in 1:100){
    vect[x] <- x*P_0(x)
  }
  k_10 <- exp(-t)*sum(vect)
  k_00 <- ((1-p)/p)*k_10
  P_1 = function(x){
    prod2 = rep.int(0,x+1)
    for(i in 1:(x+1)){
      prod2[i] <- 1+i+p*R
    }
    prod2final <- 1
    for(i in 1:(x+1)){
      prod2final <- prod2final*prod2[i]
    }
    ((p*R)^x)*(x+2)/prod2final
  }
  vect2 <- rep.int(0,100)
  for(x in 1:100){
    vect2[x] <- x*P_1(x)
  }
  k_11 <- exp(-2*t)*sum(vect2)
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

FINAL15 <- rep.int(0,100)
for(i in 1:100){
  FINAL15[i] <- final(i/100)
}

plot(FINAL3, col="red", type="l",ylim=c(0,1))
par(new=TRUE)
plot(FINAL2, col="blue", type="l", ylim=c(0,1))
par(new=TRUE)
plot(FINAL15, type="l", ylim=c(0,1))

