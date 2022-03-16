install.packages("igraph")
library("igraph")

n = 1000 #set population size
G <- sample_pa(n, power=1, m=1, directed=F) #random undirected graph G with n nodes
V(G)$id <- 1:vcount(G) #labelling each vertex for future reference
Nei <- rep.int(0,n) #empty vector for neighbours
for(i in 1:n){
  Nei[i] <- length(neighbors(G, i))
} #vector of nodes' degree in G

b <- 0.2 #beta
m <- 0.4 #gamma
#beta, gamma irrelevant in the end as we plot the ratio R1/R0 and both are multiplied by beta/gamma

R <- rep.int(0,n)
for(l in 1:n){
  x <- 0
  for(i in neighbors(G, l)){
    x <- x + (Nei[i] - 1) 
  }
  R[l] <- x/Nei[l] #k_nn(l)
}
R0 <- (b/(m*n))*sum(R) #calculated R0 as defined in essay

RND <- function(r){ #will calculate R for strategy RND
  #now need a new induced CT graph
  CT <- sample(1:n, r*n, replace=F) #random sampling of r*n individuals
  G2 <- induced_subgraph(G, CT) #induce subgraph
  Isolated = which(degree(G2)==0) 
  G1 = delete.vertices(G2, Isolated) #exclude those with no app-using neighbours
  P <- rep.int(0,n) #empty vector for neighbor-with-app numbers, k_i'
  #k_i' = 0 for non-app users
  for(i in V(G1)$id){ #app users
    t <- which(V(G1)$id==i) #find corresponding vertex in CT network
    P[i] <- length(neighbors(G1, t)) #number of app-using neighbors
  } 
  #next ten lines are the same calculation for each strategy, using k_i=Nei[i], k_i' = P[i]
  RR <- rep.int(0,n)
  for(l in 1:n){
    x <- 0
    for(i in neighbors(G, l)){
      x <- x + P[i]*(Nei[i] - 1)/Nei[i]
    }
    RR[l] <- x/Nei[l] 
  }
  R1 <- R0 - (b/(m*n))*sum(RR) #new reproduction number as defined in essay
  return(R1/R0) #new reproduction number/R0
}

vec <- seq(0,1,by=0.01)
yvec <- rep.int(0,101)
for(i in vec){
  yvec[100*i+1] <- RND(i) #calculate R/R0 for r=0,0.01,...,1
}

#now calculate R for strategy DEG
deg <- degree(G)
order <- order(deg, decreasing=T) #ordering vertices in descending order of degree

DEG <- function(r){ #will calculate R for strategy DEG
  CT2DEG <- order[1:(r*n)] #selected r*n nodes with greatest degree
  GDEG <- induced_subgraph(G, CT2DEG) #inducing this subgraph
  IsolatedDeg = which(degree(GDEG)==0) #exclude those with no app-using neighbours
  G1DEG = delete.vertices(GDEG, IsolatedDeg) 
  P <- rep.int(0,n) #empty vector for neighbor-with-app numbers, k_i'
  #k_i' = 0 for non-app users
  for(i in V(G1DEG)$id){ #app-users
    t <- which(V(G1DEG)$id==i) #find corresponding vertex in CT network
    P[i] <- length(neighbors(G1DEG, t)) #number of app-using neighbours
  }
  #same calculation as for previous strategy
  RR <- rep.int(0,n) 
  for(l in 1:n){ 
    x <- 0
    for(i in neighbors(G, l)){ 
      x <- x + P[i]*(Nei[i] - 1)/Nei[i]
    }
    RR[l] <- x/Nei[l]
  }
  R1 <- R0 - (b/(m*n))*sum(RR)
  return(R1/R0) #new reproduction number/R0
}

vecd <- seq(0,1,by=0.01)
yvecDEG <- rep.int(0,101)
for(i in vecd){
  yvecDEG[100*i+1] <- DEG(i) #calculate R/R0 for r=0,0.01,...,1
}

#now calculate R for strategy DEG
DTI <- function(r){
  V <- rep.int(1,r*n) #empty vector, will be filled with node number of each app-user
  V[1] <- order[1] #first CT download is node with greatest degree
  for(t in 1:(r*n-1)){ #each time step until r*n people have the app
    S <- NULL #starting with an empty vector for set of people who have a neighbor with the app
    for(i in V[1:t]){ #set of i people who have the app so far
      S <- union(S,neighbors(G,i))
    } #now S is a vector of all the vertices with a neighbor with the app
    S <- S[!S %in% V[1:t]] #removing everyone who already has the app
    S1 <- rep.int(0,length(S)) #this will be how many neighbors a node has with the app
    for(i in 1:length(S)){
      K <- intersect(V[1:t],neighbors(G,S[i])) #i's app-using neighbours
      S1[i] <- length(K) #number of neighbours with app
    }
    SP <- rep.int(0,length(S)) #empty vector to be probability of selection for each vertex
    for(i in 1:length(S)){
      SP[i] <- S1[i]/sum(S1) #probability as defined in essay
    }
    V[t+1] <- sample(S,1,prob=SP) #random sample with probability SP
  }
  GDTI <- induced_subgraph(G, V) #induced subgraph
  IsolatedDTI = which(degree(GDTI)==0) #exclude those with no app-using neighbours
  G1DTI = delete.vertices(GDTI, IsolatedDTI)
  P <- rep.int(0,n) #empty vector for neighbor-with-app numbers, k_i'
  #k_i' = 0 for non-app users
  for(i in V(G1DTI)$id){ #app-users
    t <- which(V(G1DTI)$id==i) #find corresponding vertex in CT network
    P[i] <- length(neighbors(G1DTI, t)) #number of app-using neighbours
  }
  #same calculation as previous strategies
  RR <- rep.int(0,n) 
  for(l in 1:n){ 
    x <- 0
    for(i in neighbors(G, l)){ 
      x <- x + P[i]*(Nei[i] - 1)/Nei[i]
    }
    RR[l] <- x/Nei[l]
  }
  R1 <- R0 - (b/(m*n))*sum(RR)
  return(R1/R0) #new reproduction number/R0
}

yvecDTI <- rep.int(0,99)
vecdti <- seq(0.01,0.99,by=0.01)
for(i in vecdti){
  yvecDTI[100*i] <- DTI(i) #calculate R/R0 for r=0.01,...,0.99
} 
yvecDTI2 <- c(1,yvecDTI,0) #add on values for r=0, r=1

pdf("ratios.pdf", width=6, height=6)
lo <- loess(yvec~seq(0,1,by=0.01)) #smoother for the random results 
plot(predict(lo),lwd=2,type="l",xlim=c(0,100),ylim=c(0,1),xlab="App uptake (%)",ylab="R/R0")
loDEG <- loess(yvecDEG~seq(0,1,by=0.01)) #smoother for the random results 
par(new=TRUE)
plot(predict(loDEG),col='blue',lwd=2,type="l",xlim=c(0,100),ylim=c(0,1),xlab="App uptake (%)",ylab="R/R0")
loDTI2 <- loess(yvecDTI2~seq(0,1,by=0.01)) #smoother for the random results 
par(new=TRUE)
plot(predict(loDTI2),col='red',lwd=2,type="l",xlim=c(0,100),ylim=c(0,1),xlab="App uptake (%)",ylab="R/R0")
legend(0.75,0.25, legend=c("RND", "DEG", "DTI"),col=c("black", "blue", "red"), lty=1:1:1, cex=0.8)
dev.off()

