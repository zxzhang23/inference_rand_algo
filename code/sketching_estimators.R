### generate sketching estimators

library(mvtnorm)
library(phangorn)
library(ggplot2)
library(tidyr)

### i.i.d. sketching
sampleDist<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/3,1/3,1/3))*sqrt(3/2)
}

sampleDist2<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/6,2/3,1/6))*sqrt(3)
}

sampleDistsparse<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/20,9/10,1/20))*sqrt(10)
}

# different types of S have different kurtosis
Esticoef_iid<-function(m,n,c,X,y,type=6,partial=0){
  if(type==1){S<-matrix(sampleDist(m*n),m)/sqrt(m)}  
  if(type==2){S<-matrix(runif(m*n,-sqrt(3),sqrt(3)),m)/sqrt(m)}
  if(type==3){S<-matrix(sampleDist2(m*n),m)/sqrt(m)}
  if(type==4){S<-matrix(rnorm(m*n),m)/sqrt(m)}
  if(type==5){S<-matrix(rt(m*n,10),m,n)/sqrt(10/8)/sqrt(m)}
  if(type==6){S<-matrix(sampleDistsparse(m*n),m)/sqrt(m)}
  
  P<-S%*%X
  M<-solve(t(P)%*%P)
  if(partial==0){w<-t(P)%*%(S%*%y);g<-M%*%w}
  else{g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g)))
}

### uniform orthogonal sketching
Generate_Haar<-function(m,n){
  O<-matrix(rnorm(m*n),m,n)
  S<-t(svd(O)$v)*sqrt(n/m)
  return(S)
}

Esticoef_Haar<-function(m,n,c,X,y,partial=0){
  S_haar<-Generate_Haar(m,n)
  P<-S_haar%*%X
  M<-solve(t(P)%*%P)
  if(partial==0){ w<-t(P)%*%(S_haar%*%y); g<-M%*%w }
  else{g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g)))
}

### Hadamard sketching
padding<-function(X,y){
  m<-nrow(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
    pady<-append(y,rep(0,2^m1-m))
  }
  else{padX<-X;pady<-y}
  return(list(padX=padX,pady=pady))
}

Esticoef_SRHT<-function(m,c,X,y,partial=0){
  p<-ncol(X)
  pad<-padding(X,y)
  X1<-pad$padX;y1<-pad$pady
  n1<-nrow(X1)
  gamma<-m/n1
  SXy<-apply(sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))*cbind(X1,y1),2,fhm)[which(rbinom(n1,1,gamma)!=0),]/sqrt(m)
  if(partial==0){g<-qr.solve(SXy[,1:p],SXy[,p+1])}
  else{SX<-SXy[,1:p];M<-solve(t(SX)%*%SX);g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g),r3=SXy))
}
