### Figure 7

library(mvtnorm)
library(phangorn)
library(ggplot2)
library(tidyr)
#setwd('C:/Users/wld301912/Desktop/test111111')
source("sketching_methods.R")

##Estimating for the largest singular value
n=2^11;p=15;k=1
grid_m=seq(800,1800,200)
set.seed(15)
#Case 1
d<-1/(1:p)
D<-diag(d)
O1<-matrix(rnorm(n*p),n,p)
W<-svd(O1)$u
O2<-matrix(rnorm(p*p),p,p)
U<-qr.Q(qr(O2))
X<-W %*% D %*% t(U)
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

Sigma<-svd(X)$d
v<-sort(Sigma,decreasing = TRUE)
sigmaj<-v[k]
lambdaj<-sigmaj^2

####################
##theoretical asymptotic variances for Hadamard sketching
th_v_hadam<-3*(1-grid_m/n)*lambdaj^2 

####################
###theoretical asymptotic variances for i.i.d. sketching estimators and Countsketch
th_v_iid<-rep(2*lambdaj^2,length(grid_m))
th_v_countsketch<-th_v_iid

####################
###theoretical asymptotic variances for Haar sketching estimators
th_v_haar<-2*(1-grid_m/n)*lambdaj^2 

#####################
### simulations
# sketch-and-solve for uniform sampling, i.i.d., and Hadamard
sim=500
#prob<-rep(1/n,n)
#r_subsample<-matrix(0,sim,length(grid_m))

#for(i in 1:sim){
#  for(j in 1:length(grid_m)){
#    ind<-sample(n,grid_m[j],replace=TRUE,prob=prob)
#    SX<-X[ind, ];Sy<-y[ind]
#    D<-(grid_m[j]*pi[ind])^(-1)
#    SXy<-(t(SX))%*%(D*Sy)
#    subls<-solve(t(SX)%*%(D*SX))%*%SXy
#    r_subsample[i,j]<-sum(c*subls)}
#}

r_iid<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
     r_iid[i,j]<-Estisingval_iid(k,grid_m[j],n,X,type=1)}
}


r_srht<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_srht[i,j]<-Estisingval_SRHT(k,grid_m[j],n,X)}
}

#r_haar<-matrix(0,sim,length(grid_m))
#for(i in 1:sim){
#  for(j in 1:length(grid_m)){
#    r_haar[i,j]<-Estisingval_Haar(k,grid_m[j],n,X)}
#}

r_countsketch<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_countsketch[i,j]<-Estisingval_CountSketch(k,grid_m[j],n,X)}
}

#################################################
### display the variances
sketch_size<-as.vector(t(matrix(rep(grid_m,sim),length(grid_m))))
df_var1<-data.frame(sketch_size,#uniform_sampling=as.vector(r_subsample)*sqrt(sketch_size),
                    iid=as.vector(r_iid)*sqrt(sketch_size),hadamard=as.vector(r_srht)*sqrt(sketch_size),countsketch=as.vector(r_countsketch)*sqrt(sketch_size))

v1<-aggregate(.~sketch_size,df_var1, FUN = stats::var)

v1$hadamard_theory<-th_v_hadam
v1$iid_theory<-th_v_iid
#v1$haar_theory<-th_v_haar
v1$countsketch_theory<-th_v_countsketch
res_vr1<-gather(v1,type,var,-sketch_size)
#res_vr1$var<-log(res_vr1$var,10)


p1f<-ggplot(res_vr1,aes(x=sketch_size/n,var,group=type,color=type,linetype=type,pointtype=type))+
  geom_point(aes(shape=type,color=type),size=1.8)+geom_line(aes(size=type))+
  scale_shape_manual(values=c(15,3,16,4,17,5))+
  scale_size_manual(values=c(0.8,0.8,0.8,0.8,0.8,0.8))+
  scale_linetype_manual(values = c(1,5,4,3,2,6))+
  scale_color_manual(values=c('#FF0000','#000000','#0066CC','green','#FFAA00','#BBCCFF'))+
  labs(title = "SVD sketching", x = 'm/n')+ylab(expression("var"))+
  theme(plot.title = element_text(hjust = 0.5),legend.box.background = element_rect(color="black", size=1),
        legend.text=element_text(size=10),axis.text.x = element_text(size = 12), axis.title.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(color = "grey20", size = 15))

#file_path<-"C:/Users/wld301912/Desktop/test111111/plots/var.pdf"
file_path<-"var.pdf"
#save(eps_var,file=file_path)
ggsave(file_path, width=4.0, height = 3.0)
