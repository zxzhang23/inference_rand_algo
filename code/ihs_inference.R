library(mvtnorm)
library(phangorn)

#################################################################################################
itera_ske<-function(m,c,X,y,ite){
  Xty<-t(X)%*%y
  beta<-matrix(0,p,ite+1)
  skedat<-matrix(0,m,p*ite)
  for(i in 2:(ite+1)){
    gamma=m/n
    Xbeta<-X%*%beta[,i-1];XtXbeta<-t(X)%*%Xbeta
    S<-matrix(rnorm(m*n),m)/sqrt(m)
    SX<-S%*%X
    M<-solve(t(SX)%*%SX)
    beta[,i]<-M%*%(Xty-XtXbeta)+beta[,i-1]
    skedat[,(1+(i-2)*p):((i-1)*p)]<-SX
  }
  return(list(r1=beta[,2:(ite+1)],r2=as.numeric(t(c)%*%beta[,2:(ite+1)]),r3=skedat))
}



######do not refresh the iteration matrices
itera_ske2<-function(m,c,X,y,ite){
  Xty<-t(X)%*%y
  beta<-matrix(0,p,ite+1)
  skedat<-matrix(0,m,p*ite)
  S<-matrix(rnorm(m*n),m)/sqrt(m)
  SX<-S%*%X
  M<-solve(t(SX)%*%SX)
  
  for(i in 2:(ite+1)){
    gamma=m/n
    Xbeta<-X%*%beta[,i-1];XtXbeta<-t(X)%*%Xbeta
    beta[,i]<-M%*%(Xty-XtXbeta)+beta[,i-1]
    skedat[,(1+(i-2)*p):((i-1)*p)]<-SX
  }
  return(list(r1=beta[,2:(ite+1)],r2=as.numeric(t(c)%*%beta[,2:(ite+1)]),r3=skedat))
}


GOE<-function(n){
  Z=matrix(rnorm(n*n),n,n);
  GOE<-(Z+t(Z))/sqrt(2)
  return(GOE)
}

############################################################################################
############### several inference methods  #################################################
###  pivotal method
pivo_ite<-function(c,sX,X,y,itebeta,alpha=0.1){
  Xty<-t(X)%*%y
  m<-nrow(sX)
  ite<-length(itebeta)
  svdsX<-svd(sX)
  U<-svdsX$u;D<-svdsX$d;V<-svdsX$v
  A<-V%*%diag(D^(-1))
  Uty<-diag(D^(-1))%*%t(V)%*%Xty
  l=matrix(0,500,ite)
  for(i in 1:500){
    GUty = Uty
    for(j in 1:ite){
      G<-GOE(p);
      GUty<-G%*%GUty;
      l[i,j]<-sum(c*(A%*%GUty))/m^(j/2)
    }
  }
  lb<-itebeta-apply(l,2,quantile,c(1-alpha/2))
  rb<-itebeta-apply(l,2,quantile,c(alpha/2))
  
  return(list(lb=lb,rb=rb))
}

###pivotal inference without refreshing sketching 
pivo_ite2<-function(c,sX,X,y,itebeta,alpha=0.1){
  Xty<-t(X)%*%y
  m<-nrow(sX)
  ite<-length(itebeta)
  svdsX<-svd(sX)
  U<-svdsX$u;D<-svdsX$d;V<-svdsX$v
  A=-V%*%diag(D^(-1))   #including the negative sign here!
  Uty<-diag(D^(-1))%*%t(V)%*%Xty
  l=matrix(0,500,ite)
  for(i in 1:500){
    GUty = Uty
    G<-GOE(p)
    for(j in 1:ite){
      GUty<-G%*%GUty;
      l[i,j]<-sum(c*(A%*%GUty))/m^(j/2)
    }
  }
  lb<-itebeta-apply(l,2,quantile,c(1-alpha/2))
  rb<-itebeta-apply(l,2,quantile,c(alpha/2))
  
  return(list(lb=lb,rb=rb))
}


### sub_rand, multi-run plug-in and, multi-run aggregation

tau<-function(m,j){m^(j/2)}
sub_multi_ite<-function(b,m,c,X,y,itebeta,K,alpha=0.05){
  ite<-length(itebeta)
  subske=matrix(0,K,ite) 
  for(i in 1:K){
    itbe<-itera_ske2(b,c,X,y,ite)$r2 ##no refresh
    subske[i,]<-itbe  
  }
  taum<-tau(m,seq(1,ite,1))
  taub<-tau(b,seq(1,ite,1))
  
  sub_lb<-itebeta-apply(t((taub/(taum-taub))*t(subske-rep(1,K)%*%t(itebeta))),2,quantile,c(1-alpha/2))
  sub_rb<-itebeta-apply(t((taub/(taum-taub))*t(subske-rep(1,K)%*%t(itebeta))),2,quantile,c(alpha/2))
  
  
  multi_v<-apply(subske,2,var)
  multi_center<-apply(subske,2,mean)
  multi_lb<-multi_center-qnorm(1-alpha/2,sd=sqrt(multi_v/K))
  multi_rb<-multi_center-qnorm(alpha/2,sd=sqrt(multi_v/K))
  
  return(list(sub_lb=sub_lb,sub_rb=sub_rb,multi_lb=multi_lb,multi_rb=multi_rb))
}



########################################################
boot_ite<-function(m,c,X,y,skedat,itebeta,K,alpha=0.1){
  Xty<-t(X)%*%y
  ite<-length(itebeta)
  r_bootstrap<-matrix(0,K,ite)
  beta=boot_beta<-matrix(0,p,ite+1)
  for(j in 2:(ite+1)){
    sX<-skedat[,(1+(j-2)*p):((j-1)*p)]
    M<-solve(t(sX)%*%sX)
    Xbeta<-X%*%beta[,j-1];XtXbeta<-t(X)%*%Xbeta
    beta[,j]<-M%*%(Xty-XtXbeta)+beta[,j-1]
    for(i in 1:K){
      smp<-sample(m,replace=TRUE)
      A<-sX[smp,]
      M<-solve(t(A)%*%A)
      Xbeta<-X%*%beta[,j-1];XtXbeta<-t(X)%*%Xbeta
      boot_beta[,j]<-M%*%(Xty-XtXbeta)+beta[,j-1]
      r_bootstrap[i,j-1]<-sum(c*boot_beta[,j])
    }
  }
  
  lb<-itebeta-apply((r_bootstrap-rep(1,K)%*%t(itebeta)),2,quantile,c(1-alpha/2))
  rb<-itebeta-apply((r_bootstrap-rep(1,K)%*%t(itebeta)),2,quantile,c(alpha/2))
  
  return(list(lb=lb,rb=rb))
}




