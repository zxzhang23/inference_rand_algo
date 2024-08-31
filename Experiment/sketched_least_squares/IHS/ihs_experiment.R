
###################################################################
# "source" the file ihs_inference.R first
compare_ite<-function(c,m,X,y,ite,b,K,sim,alpha=0.1){
  pivo_conf=sub_conf=boot_conf=multi_conf=matrix(0,sim,2*ite)
  
  for(i in 1:sim){
    res<-itera_ske2(m,c,X,y,ite)
    skedat<-res$r3  #used in bootstrap
    
    pivo_int<-pivo_ite2(c,skedat[,1:p],X,y,res$r2,alpha=0.1)
    pivo_conf[i,seq(1,2*ite,2)]<-pivo_int$lb
    pivo_conf[i,seq(2,2*ite,2)]<-pivo_int$rb
    
    
    res_two<-sub_multi_ite(b,m,c,X,y,res$r2,K,alpha=0.1)
    sub_conf[i,seq(1,2*ite,2)]<-res_two$sub_lb
    sub_conf[i,seq(2,2*ite,2)]<-res_two$sub_rb
    
    multi_conf[i,seq(1,2*ite,2)]<-res_two$multi_lb
    multi_conf[i,seq(2,2*ite,2)]<-res_two$multi_rb
    
    boot_int<-boot_ite(m,c,X,y,skedat,res$r2,K,alpha)
    boot_conf[i,seq(1,2*ite,2)]<-boot_int$lb
    boot_conf[i,seq(2,2*ite,2)]<-boot_int$rb
  }

  return(list(pivo_conf=pivo_conf,sub_conf=sub_conf,boot_conf=boot_conf,multi_conf=multi_conf))
}



p=10;n=5000
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
ls<-qr.solve(X,y)
sprintf("%.20f",ls[1])

set.seed(NULL)
c<-c(1,rep(0,p-1))

res<-compare_ite(c,1000,X,y,10,500,50,100,alpha=0.1)


pivo_conf<-res$pivo_conf
sub_conf<-res$sub_conf
boot_conf<-res$boot_conf
multi_conf<-res$multi_conf
conf<-cbind(pivo_conf,sub_conf,boot_conf,multi_conf)


decimal_places <- 15
options(digits = decimal_places)

write.csv(conf,paste("0p10n5000m1000b100to600K50confiteiid_",index,".csv",sep=""))
