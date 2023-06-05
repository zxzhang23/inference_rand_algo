library(mvtnorm)
library(phangorn)


##################################################
###inference from SRHT sketching estimators

# parameters：
# n: sample size n, p: dimension 
# m: sketching size 
# c: coefficient vector c, e.g., c= (1,0,0,\cdots,0)
# K: number of repeated experiments (and also bootstrap samples)
# b: subsketching size


### if X and y are not of size to the power of two, padding them with zeros
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


# Input: sketch size m, coefficient vector c, data X, y
# Output: sketched least square solutions r1, linear combination of r1, sketched data S*X and S*y
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

###################################################
##Pivotal method 
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

###############################################################
###Subsketching
###Input: b is subskeching size, m is sketching size, X and y are used to generate repeated samples of smaller size b
sub_int_hadamard<-function(b,m,c,X,y,sX,sy,K,alpha=0.05){
  p<-ncol(X)
  pad<-padding(X,y)
  X1<-pad$padX;y1<-pad$pady
  n1<-nrow(X1)
  gamma=b/n1
  subske=0
  for(i in 1:K){
    SXy<-Esticoef_SRHT(b,c,X,y,partial=0)$r3
    g<-qr.solve(SXy[,1:p],SXy[,p+1])
    subske[i]<-g[1]
  }
  
  taum<-sqrt((m-p)*(n1-p)/(n1-m))
  taub<-sqrt((b-p)*(n1-p)/(n1-b))
  
  lssk<-solve(t(sX)%*%sX)%*%(t(sX)%*%sy)
  center<-sum(c*lssk)
  
  lb<-quantile(taub*(subske-center),alpha/2)/(taum-taub)
  rb<-quantile(taub*(subske-center),(1-alpha/2))/(taum-taub)
  
  return(list(conf=c(center-rb,center-lb)))
}


##########################################
###bootstrap
boot_int<-function(c,sX,sy,K,alpha){
  m<-nrow(sX);p<-ncol(sX)
  r_bootstrap<-matrix(0,p,K)
  for(i in 1:K){
    smp<-sample(m,replace=TRUE)
    A<-sX[smp,];b<-sy[smp]
    r_bootstrap[,i]<-qr.solve(A,b)
  }
  bt<-t(c)%*%r_bootstrap
  
  lssk<-solve(t(sX)%*%sX)%*%(t(sX)%*%sy)
  center<-sum(c*lssk)
  
  
  rb<-quantile(bt-center,1-alpha/2)
  lb<-quantile(bt-center,alpha/2)
  return(list(conf=c(center-rb,center-lb)))
}

##################################
###plug-in method
plug_int<-function(c,m,X,y,sX,sy,K,alpha){
  p<-ncol(X);n<-nrow(X)
  plugske=0
  for(i in 1:K){
    SXy<-Esticoef_SRHT(m,c,X,y,partial=0)$r3
    g<-qr.solve(SXy[,1:p],SXy[,p+1])
    plugske[i]<-g[1]
  }
  v<-var(plugske)
  rb=qnorm(1-alpha/2,sd=sqrt(v))
  lb=qnorm(alpha/2,sd=sqrt(v))
  
  lssk<-solve(t(sX)%*%sX)%*%(t(sX)%*%sy)
  center<-sum(c*lssk)
  return(list(conf=c(center-rb,center-lb)))
}

######################################
### multi-run
multi_run<-function(c,m,X,y,sX,sy,K,alpha){
  p<-ncol(X);n<-nrow(X)
  plugske=0
  for(i in 1:K){
    SXy<-Esticoef_SRHT(m,c,X,y,partial=0)$r3
    g<-qr.solve(SXy[,1:p],SXy[,p+1])
    plugske[i]<-g[1]
  }
  v<-var(plugske)
  rb=qnorm(1-alpha/2,sd=sqrt(v/K))
  lb=qnorm(alpha/2,sd=sqrt(v/K))
  
  center<-mean(plugske)
  return(list(conf=c(center-rb,center-lb)))
}



###############################################################
### compare the performence of the above five methods
compare_methods_hadamard<-function(c,X,y,grid_m,b1,K,sim,alpha=0.1){
  p<-ncol(X);n<-nrow(X)
  pivo_conf=sub_conf=boot_conf=plug_conf=multi_run_conf=matrix(0,sim,2*length(grid_m))
  pivo_time=sub_time=boot_time=plug_time=multi_run_time=matrix(0,sim,length(grid_m))
  for(i in 1:sim){
    for(j in 1:length(grid_m)){
      SD<-Esticoef_SRHT(grid_m[j],c,X,y,partial=0)$r3
      SX<-SD[,1:p];Sy<-SD[,p+1]
      
      
      ##be careful of the sketching type in pivo_iid_growp,sub_int,and skeske_int
      Xy<-padding(X,y)
      X1<-Xy$padX;y1<-Xy$pady
      n1=nrow(X1);
      
      
      st1<-Sys.time()
      pivo_conf[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n1,SX,Sy,partial=0,alpha)$conf
      ed1<-Sys.time()
      pivo_time[i,j] <- difftime(ed1, st1, units = "secs")
      
      
      
      st2<-Sys.time()
      sub_conf[i,(2*j-1):(2*j)]<-sub_int_hadamard(b1,grid_m[j],c,X,y,SX,Sy,K,alpha)$conf
      ed2<-Sys.time()
      sub_time[i,j] <- difftime(ed2, st2, units = "secs")
      
      st3<-Sys.time()
      boot_conf[i,(2*j-1):(2*j)]<-boot_int(c,SX,Sy,K,alpha)$conf
      ed3<-Sys.time()
      boot_time[i,j]<-difftime(ed3, st3, units = "secs") 
      
      st4<-Sys.time()
      plug_conf[i,(2*j-1):(2*j)]<-plug_int(c,grid_m[j],X,y,SX,Sy,K,alpha)$conf
      ed4<-Sys.time()
      plug_time[i,j]<-difftime(ed4, st4, units = "secs")
      
      st5<-Sys.time()
      multi_run_conf[i,(2*j-1):(2*j)]<-multi_run(c,grid_m[j],X,y,SX,Sy,K,alpha)$conf
      ed5<-Sys.time()
      multi_run_time[i,j]<-difftime(ed5, st5, units = "secs")
    }
  }
  return(list(pivo_conf=pivo_conf,sub_conf=sub_conf,boot_conf=boot_conf,plug_conf=plug_conf,multi_run_conf=multi_run_conf,
              pivo_time=pivo_time,sub_time=sub_time,boot_time=boot_time,plug_time=plug_time,multi_run_time=multi_run_time))
}







