
### pivotal inference from sketching estimators of different types,
### including i.i.d., uniform orthogonal, and Hadamard
pivo_iid<-function(c,n,sX,sy,ls_sk,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  inv_sX<-solve((t(sX)%*%sX))%*%t(sX)
  H_XS<-sX%*%inv_sX
  sep<-sy-H_XS%*%sy
  if(partial==0){
    va<-sum((t(c)%*%inv_sX)^2)*sum((t(sy)%*%(diag(1,nrow(H_XS))-H_XS))^2)
    est_v<-(gamma/((gamma-xi)))*va
  }
  else{
    est_v<-((gamma-xi)/gamma)*((sum((sX%*%ls_sk)^2)*sum((t(c)%*%inv_sX)^2)+(sum(c*ls_sk))^2))
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(lb,rb))
}



pivo_haar<-function(c,n,sX,sy,ls_sk,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  inv_sX<-solve((t(sX)%*%sX))%*%t(sX)
  H_XS<-sX%*%inv_sX
  if(partial==0){
    va<-sum((t(c)%*%inv_sX)^2)*sum((t(sy)%*%(diag(1,nrow(H_XS))-H_XS))^2)
    est_v<-(gamma*(1-gamma)/((gamma-xi)*(1-xi)))*va
  }
  else{
    est_v<-((1-gamma)*(gamma-xi)/(gamma*(1-xi)^3))*((sum((sX%*%ls_sk)^2)*sum((t(c)%*%inv_sX)^2)+(sum(c*ls_sk))^2))
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(lb,rb))
}


pivo_hadamard<-function(c,n,sX,sy,ls_sk,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  inv_sX<-solve((t(sX)%*%sX))%*%t(sX)
  H_XS<-sX%*%inv_sX
  if(partial==0){
    est_v<-((gamma*(1-gamma)/((gamma-xi)*(1-xi)))*sum((t(c)%*%inv_sX)^2)*sum((t(sy)%*%(diag(1,nrow(H_XS))-H_XS))^2))
  }
  else{
    est_v<-((1-gamma)*(gamma-xi)/(gamma*(1-xi)^3))*((sum((sX%*%ls_sk)^2)*sum((t(c)%*%inv_sX)^2)+2*(sum(c*ls_sk))^2))
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(lb,rb))
}
