args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  index <- as.numeric(args[1])  #folder number
  set.seed(index)
} else {
  stop()
}

library(phangorn) #fhm function
library(memuse)  #memory report

### if the sample size is not to the power of two, padding X,y with zeros makes them comfortable with SRHT
padding<-function(X,y){
  m<-nrow(X);p<-ncol(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
    pady<-append(y,rep(0,2^m1-m))
  }
  else{padX<-X;pady<-y}
  return(list(padX=padX,pady=pady))
}

SRHT<-function(select,D,a){
  Da<-D*a
  Sa<-fhm(Da)[select]  
  return(Sa)
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

set.seed(123)
n=1000000;p=2000
set.seed(123)
X<-matrix(rnorm(n*p),n,p)
beta<-seq(1,p)/p
y<-X%*%beta+as.vector(rnorm(n,sd=1))
set.seed(NULL)

grid_m=300000
grid_b=4000
grid_K=50


print(paste('The amount of total ram used by the current R process is', Sys.procmem()$size, 'and the maximum memory usage is', Sys.procmem()$peak))

mem_bf<-Sys.procmem()$size
n<-nrow(X);p<-ncol(X)
st<-Sys.time()
ls<-solve(qr(X,LAPACK=TRUE),y)
ed<-Sys.time()
mem_af<-Sys.procmem()$size
lstim <- difftime(ed, st, units = "secs")
paste('ram used in computing full LS',mem_af-mem_bf)
print(paste('running time for full least squares:',difftime(ed,st,units='secs'),'secs'))
print(paste('The amount of total ram used by the current R process is', Sys.procmem()$size, 'and the maximum memory usage is', Sys.procmem()$peak))
cat('\n\n')
print(ls[p])



pad<-padding(X,y)
X1<-pad$padX;y1<-pad$pady
rm(X);rm(y)
n1<-nrow(X1)
sim=4   #number of replication in each parallel job
alpha=0.1  #confidence level
pivo_conf=sub_conf=boot_conf=plug_conf=multi_run_conf=multi_run_conf1=matrix(0,sim,2*length(grid_m))
pivo_time=sub_time=boot_time=plug_time=multi_run_time=matrix(0,sim,length(grid_m))
sub_ske_all=matrix(0,sim,grid_K)
preli_time = matrix(0,sim,length(grid_m))
preli_est=0


c<-rep(0,p) 
c[p]<-1 #coefficient vector, the c[p] is one since we aim to inference on p-th coord of ls solution


for(i in 1:sim){
  for(j in 1:length(grid_m)){
    m<-grid_m[j]
    b<-grid_b[j]
    K<-grid_K[j]
    st_sxm<-Sys.time()
    gamma<-m/n1
    select0<-which(rbinom(n1,1,gamma)!=0)
    D0<-sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))  
    SXm<-matrix(0,length(select0),p)
    for(k in 1:p){
      SXm[,k]<-SRHT(select0,D0,X1[,k])
    }
    
    Sym<-SRHT(select0,D0,y1)
    ed_sxm<-Sys.time()
    
    print(paste('sim',i, 'time and ram for preliminary estimator (size m)'))
    print(paste('running time for generating sketch:',difftime(ed_sxm,st_sxm,units='secs'),'secs'))
    print(paste('the amount of total ram used by the current R process is', Sys.procmem()$size, 'and the maximum memory usage is', Sys.procmem()$peak))
    print(paste('the maximum memory usage is', Sys.procmem()$peak))
    cat('\n')
    
    print('time and ram for sketched LS part (size m)')
    st_lsm<-Sys.time()
    lsskm<-solve(qr(SXm,LAPACK=TRUE),Sym)
    ed_lsm<-Sys.time()
    print(paste('running time for solving sketched LS (size m):',difftime(ed_lsm,st_lsm,units='secs'),'secs'))
    cat('\n\n')
    
    
    preli_time[i,j] <- difftime(ed_sxm,st_sxm,units='secs')+difftime(ed_lsm, st_lsm, units = "secs")
    preli_est[i]<-sum(c*lsskm)
    
    
    st1<-Sys.time()
    pivo_conf[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n1,SXm,Sym,partial=0,alpha)$conf
    ed1<-Sys.time()
    pivo_time[i,j] <- difftime(ed1, st1, units = "secs")
    
    
    print('time and ram for subrandomization')
    mem_bf<-Sys.procmem()$size
    st2<-Sys.time()
    subske=0
    for(l in 1:K){
      gamma<-b/n1
      select<-which(rbinom(n1,1,gamma)!=0)
      D<-sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))  
      SX<-matrix(0,length(select),p)
      for(k in 1:p){
        SX[,k]<-SRHT(select,D,X1[,k])
      }
      
      Sy<-SRHT(select,D,y1)
      subske[l]<-sum(c*solve(qr(SX,LAPACK=TRUE),Sy))     
    }
    sub_ske_all[i,]<-subske
    
    
    taum<-sqrt((m-p)*(n1-p)/(n1-m))
    taub<-sqrt((b-p)*(n1-p)/(n1-b))
    center<-sum(c*lsskm)
    lb<-quantile(taub*(subske-center),alpha/2)/(taum-taub)
    rb<-quantile(taub*(subske-center),(1-alpha/2))/(taum-taub)
    sub_conf[i,(2*j-1):(2*j)]<-c(center-rb,center-lb)
    ed2<-Sys.time()
    mem_af<-Sys.procmem()$size
    paste('ram used in subrandomization',mem_af-mem_bf)
    print(paste('the maximum memory usage is', Sys.procmem()$peak))
    cat('\n\n')
    sub_time[i,j] <- difftime(ed2, st2, units = "secs")
    
    
    
    plug_v<-var(subske)
    plug_rb=qnorm(1-alpha/2,sd=sqrt(plug_v)*taub/taum)
    plug_lb=qnorm(alpha/2,sd=sqrt(plug_v)*taub/taum)
    plug_conf[i,(2*j-1):(2*j)]=c(center-plug_rb,center-plug_lb)
    
    multi_v<-var(subske)
    multi_rb=qnorm(1-alpha/2,sd=sqrt(multi_v/K))
    multi_lb=qnorm(alpha/2,sd=sqrt(multi_v/K))
    multi_center<-mean(subske)
    multi_run_conf[i,(2*j-1):(2*j)]=c(multi_center-multi_rb,multi_center-multi_lb)
  }
}

tim<-cbind(pivo_time,sub_time)
conf<-cbind(pivo_conf,sub_conf,plug_conf,multi_run_conf)



write.csv(preli_est,paste("0large_simdata_lapack_preli_est_n1e6p2e3b4e3m3e5K50rep100_",index,".csv",sep=""))
write.csv(preli_time,paste("0large_simdata_lapack_preli_time_n1e6p2e3b4e3m3e5K50rep100_",index,".csv",sep=""))
write.csv(tim,paste("0large_simdata_lapack_time_n1e6p2e3b4e3m3e5K50rep100_",index,".csv",sep=""))
write.csv(conf,paste("0large_simdata_lapack_conf_n1e6p2e3b4e3m3e5K50rep100_",index,".csv",sep=""))
write.csv(sub_ske_all,paste("0large_simdata_lapack_subskeall_n1e6p2e3b4e3m3e5K50rep100_",index,".csv",sep=""))



