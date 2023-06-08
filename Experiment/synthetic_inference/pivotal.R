library(mvtnorm)
library(phangorn)


#Case 2
p=500;n=2000
genesig<-function(t,n){
 A<-matrix(0,n,n)
 for(i in 1:n){
   for(j in 1:n){
     A[i,j]=t^(abs(i-j))
   }
 }
 return(A)
}
set.seed(15)
sigX<-2*genesig(0.5,p)
W<-rmvt(n, sigX, df = 2)
U<-qr.Q(qr(W));
V<-svd(matrix(rnorm(p*p),p))$v
ev<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
X<-U%*%ev%*%V
beta<-c(rep(1,floor(0.2*p)),0.1*rep(1,p-2*floor(0.2*p)),rep(1,floor(0.2*p)))
y<-X%*%beta+as.vector(rnorm(n,sd=0.01))
Xty<-t(X)%*%y
ls<-qr.solve(X,y)
c=c(1,rep(0,p-1))
#c=c(1,-1,rep(0,p-2))
print(sum(c*ls))
set.seed(NULL)


grid_m=seq(800,1200,200)
sim=500

pivo_conf_iid_s=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_iid(grid_m[j],n,c,X,y,type=6,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_iid_s[i,(2*j-1):(2*j)]<-pivo_iid(c,n,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_iid_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_iid(grid_m[j],n,c,X,y,type=6,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_iid_pa[i,(2*j-1):(2*j)]<-pivo_iid(c,n,SX,Sy,partial=1,alpha=0.05)$conf
  }
}


pivo_conf_haar_s=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_Haar(grid_m[j],n,c,X,y,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_haar_s[i,(2*j-1):(2*j)]<-pivo_haar(c,n,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_haar_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_Haar(grid_m[j],n,c,X,y,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_haar_pa[i,(2*j-1):(2*j)]<-pivo_haar(c,n,SX,Sy,partial=1,alpha=0.05)$conf
  }
}


Xy<-padding(X,y)
X0<-Xy$padX;y0<-Xy$pady
n0=nrow(X0)
pivo_conf_hadamard_s=matrix(0,sim,2*length(grid_m))

for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_SRHT_fast(grid_m[j],n0,c,X0,y0,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_hadamard_s[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n0,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_hadamard_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-Esticoef_SRHT_fast(grid_m[j],n0,c,X0,y0,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_hadamard_pa[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n0,SX,Sy,partial=1,alpha=0.05)$conf
  }
}

conf<-cbind(pivo_conf_iid_s,pivo_conf_haar_s,pivo_conf_hadamard_s,pivo_conf_iid_pa,pivo_conf_haar_pa,pivo_conf_hadamard_pa)

even<-seq(2,ncol(conf),2);odd<-seq(1,ncol(conf)-1,2)
right<-conf[,even];left<-conf[,odd]
accept<-(left<ls[1])&(ls[1]<right)
cov0<-apply(accept,2,sum)

cov<-matrix(cov0,ncol(conf)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio<-cov/sim
