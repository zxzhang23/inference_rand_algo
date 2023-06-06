library(mvtnorm)
library(phangorn)

##################################################
### inference from partial sketching estimators
### the main difference from the sketch-and-solve case is that $X^\top y$ is known and the partial sketching estimators should be debiased


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

## Pivotal method 
pivo_hadamard<-function(c,n,sX,sy,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  invsX<-solve((t(sX)%*%sX))
  if(partial==0){
    lssk<-invsX%*%(t(sX)%*%sy)
    center<-sum(c*lssk)
    sep<-sy-sX%*%lssk
    est_v<-(gamma*(1-gamma)/((gamma-xi)*(1-xi)))*sum(c*(invsX%*%c))*sum(sep^2)
  }
  else{
    lssk<-invsX%*%Xty
    center<-((gamma-xi)/(gamma*(1-xi)))*sum(c*lssk)
    est_v<-((1-gamma)*(gamma-xi)/(gamma*(1-xi)^3))*((sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+2*(sum(c*lssk))^2))
  }
  
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}


### Subsketching
### Note that the partial sketching estimators should be rescaled properly
sub_int_hadamard_pa<-function(b,m,c,X,y,sX,sy,K,alpha=0.05){
  p<-ncol(X)
  Xty<-t(X)%*%y
  pad<-padding(X,y)
  X1<-pad$padX;y1<-pad$pady
  n1<-nrow(X1)
  subske=0
  for(i in 1:K){
    subske[i]<-Esticoef_SRHT(b,c,X,y,partial=1)$r2
  }
  
  taum<-sqrt((m-p)*(n1-p)/(n1-m))
  taub<-sqrt((b-p)*(n1-p)/(n1-b))
  
  debiasm<-(n1*(m-p))/(m*(n1-p))
  debiasb<-(n1*(b-p))/(b*(n1-p))
  
  lssk<-debiasm*solve(t(sX)%*%sX)%*%Xty
  center<-sum(c*lssk)
  
  lb<-quantile(taub*(debiasb*subske-center),alpha/2)/(taum-taub)
  rb<-quantile(taub*(debiasb*subske-center),(1-alpha/2))/(taum-taub)
  
  return(list(conf=c(center-rb,center-lb)))
}

### bootstrap
boot_int_pa<-function(c,sX,sy,Xy,K,alpha){
  m<-nrow(sX);p<-ncol(sX)
  r_bootstrap<-matrix(0,p,K)
  for(i in 1:K){
    smp<-sample(m,replace=TRUE)
    A<-sX[smp,];
    r_bootstrap[,i]<-solve(t(A)%*%A)%*%Xy
  }
  bt<-t(c)%*%r_bootstrap
  
  lssk<-solve(t(sX)%*%sX)%*%Xy
  center<-sum(c*lssk)
  
  rb<-quantile(bt-center,1-alpha/2)
  lb<-quantile(bt-center,alpha/2)
  return(list(conf=c(center-rb,center-lb)))
}

## plug-in 
## partial sketching estimators are rescaled such that they are unbiased
plug_int_pa<-function(c,m,X,y,sX,sy,K,alpha){
  p<-ncol(X);n<-nrow(X)
  plugske=0
  for(i in 1:K){
    plugske[i]<-Esticoef_SRHT(m,c,X,y,partial=1)$r2
  }
  lssk<-solve(t(sX)%*%sX)%*%(t(X)%*%y)
  center<-sum(c*lssk)*(n*(m-p))/(m*(n-p))
  v<-var(plugske*(n*(m-p))/(m*(n-p)))
  rb=qnorm(1-alpha/2,sd=sqrt(v))
  lb=qnorm(alpha/2,sd=sqrt(v))
  
  return(list(conf=c(center-rb,center-lb)))
}




### multi-run inference
multi_run_pa<-function(c,m,X,y,sX,sy,K,alpha){
  p<-ncol(X);n<-nrow(X)
  Xty<-t(X)%*%y
  plugske=0
  for(i in 1:K){
    res<-Esticoef_SRHT(m,c,X,y,partial=1)
    plugske[i]<-res$r2
  }
  center<-mean(plugske)*(n*(m-p))/(m*(n-p))
  v<-var(plugske*(n*(m-p))/(m*(n-p)))
  rb=qnorm(1-alpha/2,sd=sqrt(v/K))
  lb=qnorm(alpha/2,sd=sqrt(v/K))
  return(list(conf=c(center-rb,center-lb)))
}



p=50;n=8000
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
Xty<-as.vector(t(X)%*%y)
ls<-solve(t(X)%*%X)%*%Xty
print(ls[1])



set.seed(NULL)
c<-c(1,rep(0,p-1))



grid_m=seq(200,2000,200)



compare_methods_hadamard_partial<-function(c,X,y,grid_m,b,K=20,sim,alpha=0.1){
  p<-ncol(X);n<-nrow(X)
  pivo_conf=sub_conf=skeske_conf=boot_conf=plug_conf=multi_run_conf=matrix(0,sim,2*length(grid_m))
  for(i in 1:sim){
    for(j in 1:length(grid_m)){
      SD<-Esticoef_SRHT(grid_m[j],c,X,y,partial=1)$r3
      SX<-SD[,1:p];Sy<-SD[,p+1]
      
      Xy<-padding(X,y)
      X1<-Xy$padX;y1<-Xy$pady
      n1=nrow(X1);
      
      pivo_conf[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n1,SX,Sy,partial=1,alpha)$conf
    
      sub_conf[i,(2*j-1):(2*j)]<-sub_int_hadamard_pa(b,grid_m[j],c,X,y,SX,Sy,K,alpha)$conf
        
      boot_conf[i,(2*j-1):(2*j)]<-boot_int_pa(c,SX,Sy,Xty,K,alpha)$conf
       
      plug_conf[i,(2*j-1):(2*j)]<-plug_int_pa(c,grid_m[j],X,y,SX,Sy,K,alpha)$conf
     
      multi_run_conf[i,(2*j-1):(2*j)]<-multi_run_pa(c,grid_m[j],X,y,SX,Sy,K,alpha)$conf
      }
  }
  return(list(pivo_conf=pivo_conf,sub_conf=sub_conf,boot_conf=boot_conf,plug_conf=plug_conf,multi_run_conf=multi_run_conf))
}


res_hadamard<-compare_methods_hadamard_partial(c,X,y,grid_m,100,K=100,sim=50,0.1)
conf<-cbind(res_hadamard$pivo_conf,res_hadamard$sub_conf,res_hadamard$boot_conf,res_hadamard$plug_conf,res_hadamard$multi_run_conf)





