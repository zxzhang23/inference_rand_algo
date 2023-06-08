
### pivotal inference from sketching estimators of different types,
### including i.i.d., uniform orthogonal, and Hadamard
### for partial sketching, assume that $X^\top y$ is known
pivo_iid<-function(c,n,sX,sy,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  invsX<-solve((t(sX)%*%sX))
  if(partial==0){
    lssk<-invsX%*%(t(sX)%*%sy)
    center<-sum(c*lssk)
    sep<-sy-sX%*%lssk
    est_v<-(gamma/((gamma-xi)))*sum(c*(invsX%*%c))*sum(sep^2)
  }
  else{
    lssk<-invsX%*%Xty
    center<-(m-p)*sum(c*lssk)/m
    est_v<-(gamma/(gamma-xi))*((sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+(sum(c*lssk))^2))
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}

pivo_haar<-function(c,n,sX,sy,partial=0,alpha){
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
    est_v<-((1-gamma)*(gamma-xi)/(gamma*(1-xi)^3))*((sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+(sum(c*lssk))^2))
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(center-rb,center-lb))
}


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
    lssk<-invsX%*%Xty    ### for partial case, we assume the knowledge of $X^\top y$
    center<-((gamma-xi)/(gamma*(1-xi)))*sum(c*lssk)
    est_v<-((1-gamma)*(gamma-xi)/(gamma*(1-xi)^3))*((sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+2*(sum(c*lssk))^2))  
    ### the asymptotic variance of Hadamard partial sketching estimators are slightly different from that of Haar case
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(center-rb,center-lb))
}
