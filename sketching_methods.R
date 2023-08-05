library(phangorn)

#iid sketching
Estisingval_iid<-function(i,m,n,X,type=2){
  if(type==1){S<-matrix(rnorm(m*n),m,n)/sqrt(m)}
  if(type==2){nbinom<-6
  p<-0.5+1/(2*sqrt(3))
  S<-matrix(rbinom(m*n,nbinom,p)-nbinom*p,m,n)/sqrt(m)}
  P<-S%*%X
  lambda<-svd(P)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmai=v[i]
  return(sigmai^2)
}

#random orthogonal sketching
Generate_Haar<-function(m,n){
  O<-matrix(rnorm(m*n),m,n)
  S<-t(svd(O)$v)*sqrt(n/m)
  return(S)
}

Estisingval_Haar<-function(i,m,n,X){
  S_haar<-Generate_Haar(m,n)
  P<-S_haar%*%X
  lambda<-svd(P)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmai=v[i]
  return(sigmai^2)
}


### Hadamard sketching
padding<-function(X){
  m<-nrow(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
  }
  else{padX<-X}
  return(padX)
}

Estisingval_SRHT<-function(i,m,n,X){
  p<-ncol(X)
  X1<-padding(X)
  n1<-nrow(X1)
  gamma<-m/n1
  SX<-apply(sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))*X1,2,fhm)[which(rbinom(n1,1,gamma)!=0),]/sqrt(m)
  lambda<-svd(SX)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmai=v[i]
  return(sigmai^2)
}

### CountSketch
Estisingval_CountSketch<-function(i,m,n,X){
  p<-ncol(X)
  SX<-matrix(0,m,p)
  X1<-sample(c(1,-1),n,replace=TRUE,prob=c(0.5,0.5))*X
  hash<-sample(m,n,replace=TRUE)
  for(j in 1:m){
    indicesj<-which(hash==j)
    if(length(indicesj) == 0){
      SX[j,]<-0 
    } else if (length(indicesj) == 1){
      SX[j,]<-X1[indicesj,]
    } else{
      SX[j,]<-colSums(X1[indicesj,])}
  }
  lambda<-svd(SX)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmai=v[i]
  return(sigmai^2)
}
