###read the functions from "inference_methods_hadamard.R" and also HGDP data
setwd("C:/Users/zhixz/Dropbox (Personal)/R/sketching/code_github")
source("inference_methods_hadamard.R")


genesig<-function(t,n){
  A<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j]=t^(abs(i-j))
    }
  }
  return(A)
}

p=500;n=8000
set.seed(15)
sigX<-2*genesig(0.5,p)
W<-rmvt(n, sigX, df = 2)
U<-qr.Q(qr(W));
V<-svd(matrix(rnorm(p*p),p))$v
ev<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
X<-U%*%ev%*%V
beta<-c(rep(1,floor(0.2*p)),0.1*rep(1,p-2*floor(0.2*p)),rep(1,floor(0.2*p)))
y<-X%*%beta+as.vector(rnorm(n,sd=0.01))
ls<-qr.solve(X,y)
c<-c(1,rep(0,p-1))
print(sum(c*ls))

grid_m<-seq(1500,4500,300)
res_hadamard<-compare_methods_hadamard(c,X,y,grid_m,1000,K=100,sim=500,alpha=0.1)
conf<-cbind(res_hadamard$pivo_conf,res_hadamard$sub_conf,res_hadamard$boot_conf,res_hadamard$plug_conf,res_hadamard$multi_run_conf)
write.csv(conf,"0p500n8000b1000K100_conf.csv")
