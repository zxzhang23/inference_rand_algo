
###read the functions from "inference_methods_hadamard.R" and also HGDP data
setwd("C:/Users/zhixz/Dropbox (Personal)/R/sketching/code_github")
source("inference_methods_hadamard.R")

url <- "https://raw.githubusercontent.com/dobriban/DPA/master/Experiments/Experiment%2011%20-%20HGDP/hgdp/chr22/hgdp.txt"
file_path <- "C:/Users/zhixz/Dropbox (Personal)/R/sketching/code_github/hgdp.txt"  

download.file(url, file_path)
hgdp<-read.table(file_path,sep=',')


####  select the first 200 features to form X, and the next feature as y
X0<-as.matrix(hgdp[,1:200])
y0<-as.matrix(hgdp[,201])
X<-scale(X0);y<-scale(y0)
n<-nrow(X);p<-ncol(X)
ls<-qr.solve(X,y)
c<-c(1,rep(0,p-1))
print(sum(c*ls))


grid_m=seq(400,800,50)
res_hadamard<-compare_methods_hadamard(c,X,y,grid_m,300,K=100,sim=5,alpha=0.1)
conf<-cbind(res_hadamard$pivo_conf,res_hadamard$sub_conf,res_hadamard$boot_conf,res_hadamard$plug_conf,res_hadamard$multi_run_conf)
write.csv(conf,"0hgdp_conf.csv")




