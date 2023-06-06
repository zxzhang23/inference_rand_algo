### Figure 7

library(mvtnorm)
library(phangorn)
library(ggplot2)
library(tidyr)

n=2^11;p=500
grid_m=seq(800,1800,200)
set.seed(15)
#Case 1
X<-matrix(rnorm(p*n),n);y<-runif(n,0,1);tXy<-t(X)%*%y
set.seed(NULL)

#Case 2
#genesig<-function(t,n){
#  A<-matrix(0,n,n)
#  for(i in 1:n){
#    for(j in 1:n){
#      A[i,j]=t^(abs(i-j))
#    }
#  }
#  return(A)
#}
#sigX<-2*genesig(0.5,p)
#W<-rmvt(n, sigX, df = 2)
#U<-qr.Q(qr(W));V<-svd(matrix(rnorm(p*p),p))$v
#ev<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
#X<-U%*%ev%*%V
#beta<-c(rep(1,floor(0.2*p)),0.1*rep(1,p-2*floor(0.2*p)),rep(1,floor(0.2*p)))
#y<-X%*%beta+as.vector(rnorm(n,sd=0.01))

c<-c(1,rep(0,p-1))
lsbeta<-qr.solve(X,y)
MSS<-sum((X%*%lsbeta)^2)
TSS<-sum(y^2)
snr<-MSS/(TSS-MSS);snr

####################
##theoretical asymptotic variances for Hadamard sketching, 
## predicted by the results for uniform orthogonal sketching
U1<-svd(X)$u
t1<-grid_m*(n-grid_m)/((n-p)*(grid_m-p))
l<-t(c)%*%Xinv
t2<-sum(l^2)*(sum(y^2)-sum((t(U1)%*%y)^2))
th_v<-t1*t2;th_v #for sketch-and-solve estimators

t3<-sum((t(U1)%*%y)^2)*sum(l^2)*(grid_m/n)^2
t4<-grid_m*(n-grid_m)*(n-p)/(grid_m-p)^3*(1+(th_m^2/(sum((t(U1)%*%y)^2)*sum(l^2))))
th_v2<-t3*t4;th_v2 #for partial sketching estimators

####################
###theoretical asymptotic variances for i.i.d. sketching estimators
kappa4 = 10
ep<- as.vector(y - X%*%lsbeta)
th_v_iid<-(kappa4 - 3)*sum((l*ep)^2)+ (grid_m/(grid_m-p))*sum(l^2)*sum(ep^2)
th_v2_iid<-(grid_m/(grid_m-p))^2*(kappa4 - 3)*sum((l*(y-ep))^2)+(grid_m/(grid_m-p))^3*(sum(l^2)*sum((y-ep)^2)+sum(c*lsbeta)^2)


#####################
### simulations
# sketch-and-solve for uniform sampling, i.i.d., and Hadamard
sim=500
prob<-rep(1/n,n)
r_subsample<-matrix(0,sim,length(grid_m))

for(i in 1:sim){
  for(j in 1:length(grid_m)){
    ind<-sample(n,grid_m[j],replace=TRUE,prob=prob)
    SX<-X[ind, ];Sy<-y[ind]
    D<-(grid_m[j]*pi[ind])^(-1)
    SXy<-(t(SX))%*%(D*Sy)
    subls<-solve(t(SX)%*%(D*SX))%*%SXy
    r_subsample[i,j]<-sum(c*subls)}
}

r_iid<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_iid[i,j]<-Esticoef_iid(grid_m[j],n,c,X,y,type=6)$r2}
}


r_srht<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_srht[i,j]<-Esticoef_SRHT(grid_m[j],c,X,y,partial=0)$r2}
}


# partial sketching

r_subsample_p<-matrix(0,sim,length(grid_m))

for(i in 1:sim){
   for(j in 1:length(grid_m)){
      ind<-sample(n,grid_m[j],replace=TRUE,prob=pi)
      SX<-X[ind, ]
      D<-(grid_m[j]*pi[ind])^(-1)
      subls<-solve(t(SX)%*%(D*SX))%*%(t(X)%*%y)
      r_subsample_p[i,j]<-sum(c*subls)}
}


r_iid_p<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
      r_iid_p[i,j]<-Esticoef_iid(grid_m[j],n,c,X,y,type=2,partial=1)$r2}
}



r_srht_p<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
      r_srht_p[i,j]<-Esticoef_SRHT(grid_m[j],c,X,y,partial=1)$r2}
}

#################################################
### display the variances
sketch_size<-as.vector(t(matrix(rep(grid_m,sim),length(grid_m))))
df_var1<-data.frame(sketch_size,uniform_sampling=as.vector(r_subsample)*sqrt(sketch_size),iid=as.vector(r_iid)*sqrt(sketch_size),
                    hadamard=as.vector(r_srht)*sqrt(sketch_size))
df_var2<-data.frame(uniform_sampling=as.vector(r_subsample_p)*sqrt(sketch_size),iid=as.vector(r_iid_p)*sqrt(sketch_size),hadamard=as.vector(r_srht_p)*sqrt(sketch_size))

v1<-aggregate(.~sketch_size,df_var1, FUN = stats::var)
v2<-aggregate(.~sketch_size,df_var2, FUN = stats::var)

v1$hadamard_theory<-th_v;v2$hadamard_theory<-th_v2
v1$iid_theory<-th_v_iid;v2$iid_theory<-th_v2_iid
res_vr1<-gather(v1,type,var, -sketch_size)
res_vr2<-gather(v2,type,var, -sketch_size)
res_vr1$var<-log(res_vr1$var,10)
res_vr2$var<-log(res_vr2$var,10)


p1f<-ggplot(res_vr1,aes(x=sketch_size/n,var,group=type,color=type,linetype=type,pointtype=type))+geom_point(aes(shape=type,color=type),size=1.8)+geom_line(aes(size=type))+
  scale_shape_manual(values=c(15,3,16,4,17))+
  scale_size_manual(values=c(0.8,0.8,0.8,0.8,0.8))+
  scale_linetype_manual(values = c(1,5,4,3,2))+
  scale_color_manual(values=c('#FF0000','#000000','#0066CC','green','#FFAA00'))+
  labs(title = "Complete sketching", x = 'm/n')+ylab(expression(log[10]("var")))+theme(plot.title = element_text(hjust = 0.5),legend.box.background = element_rect(color="black", size=1),legend.text=element_text(size=10),axis.text.x = element_text(size = 12), axis.title.x = element_text(color = "grey20", size = 15),axis.text.y = element_text(size = 12), axis.title.y = element_text(color = "grey20", size = 15))


p1p<-ggplot(res_vr2,aes(x=sketch_size/n,var,group=type,color=type,linetype=type,pointtype=type))+geom_point(aes(shape=type,color=type),size=1.8)+geom_line(aes(size=type))+
  scale_shape_manual(values=c(15,3,16,4,17))+
  scale_size_manual(values=c(0.8,0.8,0.8,0.8,0.8))+
  scale_linetype_manual(values = c(1,5,4,3,2))+
  scale_color_manual(values=c('#FF0000','#000000','#0066CC','green','#FFAA00'))+
  labs(title = "Partial sketching", x = 'm/n')+ylab(expression(log[10]("var")))+theme(plot.title = element_text(hjust = 0.5),legend.box.background = element_rect(color="black", size=0.5),legend.text=element_text(size=10),axis.text.x = element_text(size = 12), axis.title.x = element_text(color = "grey20", size = 15),axis.text.y = element_text(size = 12), axis.title.y = element_text(color = "grey20", size = 15))


combined <- p1f + p1p & theme(legend.position = "right",legend.title=element_blank())
combined + plot_layout(guides = "collect")

eps_var<-recordPlot()
file_path<-"C:/Users/zhixz/Dropbox (Personal)/R/sketching/code_github/var.eps"
save(eps_var,file=file_path)


