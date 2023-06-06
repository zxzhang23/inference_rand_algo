###############################################
## display the coverage and the length of confidence intervals

library(PropCIs)
library(ggplot2)
#############
#load the results, grid_m, and the first coordinate of the real solution
conf0= read.csv("C:/Users/zhixz/Dropbox (Personal)/R/sketching/code_randomized_algorithms/real data analysis/0hgdp_conf.csv", header =TRUE, sep = ",")
conf=conf0[,-1]
sim=nrow(conf)
grid_m=seq(400,800,50)
ls[1]<-0.03470465

#######################

even<-seq(2,ncol(conf),2);odd<-seq(1,ncol(conf)-1,2)
right<-conf[,even];left<-conf[,odd]
accept<-(left<ls[1])&(ls[1]<right)
cov0<-apply(accept,2,sum)

cov<-matrix(cov0,ncol(conf)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio0<-cov/sim
ratio<-ratio0

ci<-function(x){exactci(x,sim,0.95)$conf.int[1:2]}

int=matrix(0,ncol(conf)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    int[i,(2*j-1):(2*j)]=round(ci(cov[i,j]),digits=3)
  }
}

lower<-int[,seq(1,2*length(grid_m),2)]
upper<-int[,seq(2,2*length(grid_m),2)]

library(ggplot2)
x <- grid_m
df_full <- data.frame(x = grid_m, y1 =ratio[1,], y2=ratio[2,],y3=ratio[3,],y4=ratio[4,],y5=ratio[5,], y_upper1 =upper[1,] , y_lower1 = lower[1,],y_upper2=upper[2,],y_lower2=lower[2,],y_upper3=upper[3,],y_lower3=lower[3,],y_upper4=upper[4,],y_lower4=lower[4,],y_upper5=upper[5,],y_lower5=lower[5,])


len<-apply(right-left,2,mean)
conflen<-matrix(len,5,length(grid_m),byrow=TRUE)
len_sd<-apply(right-left,2,sd)
lensd<-matrix(len_sd,5,length(grid_m),byrow=TRUE)
low<-conflen-lensd
upp<-conflen+lensd

df_len<-data.frame(x=grid_m,y1=conflen[1,],y2=conflen[2,],y3=conflen[3,],y4=conflen[4,],y5=conflen[5,],low1=low[1,],low2=low[2,],low4=low[4,],low5=low[5,],upp1=upp[1,],upp2=upp[2,],upp4=upp[4,],upp5=upp[5,])


ggplot(df_full, aes(x = x)) +
  ylim(0.8, 1) +
  geom_point(aes(y = y1, color = "Pivotal",shape='Pivotal'), size = 2) +
  geom_line(aes(y = y1, color = "Pivotal",linetype = 'Pivotal')) +
  geom_point(aes(y = y2, color = "Subsketching",shape="Subsketching"), size = 2) +
  geom_line(aes(y = y2, color = "Subsketching",linetype = 'Subsketching')) +
  geom_point(aes(y = y3, color = "Bootstrap",shape="Bootstrap"), size = 2) +
  geom_line(aes(y = y3, color = "Bootstrap",linetype = 'Bootstrap')) +
  geom_point(aes(y = y4, color = "Plug-in",shape='Plug-in'), size = 2) +
  geom_line(aes(y = y4, color = "Plug-in",linetype='Plug-in')) +
  geom_point(aes(y = y5, color = "Multi-run",shape='Multi-run'), size = 2) +
  geom_line(aes(y = y5, color = "Multi-run", linetype = 'Multi-run')) +
  geom_errorbar(aes(ymin = lower[1,], ymax = upper[1,],color = 'Pivotal'), width =  0.05*x[1]) +
  geom_errorbar(aes(ymin = lower[2,], ymax = upper[2,],color='Subsketching'), width =  0.05*x[1]) +
  geom_errorbar(aes(ymin = lower[3,], ymax = upper[3,],color='Bootstrap'), width =  0.05*x[1]) +
  geom_errorbar(aes(ymin = lower[4,], ymax = upper[4,],color='Plug-in'), width =  0.05*x[1]) +
  geom_errorbar(aes(ymin = lower[5,], ymax = upper[5,],color='Multi-run'), width =  0.05*x[1]) +
  scale_shape_manual(values = c("Pivotal" = 1, "Subsketching" = 2, "Bootstrap" = 3, "Plug-in" = 4,"Multi-run"=5)) +
  scale_color_manual(values = c("Pivotal" = 'black', "Subsketching" = 'red', "Bootstrap" = 'blue', "Plug-in" = '355E3B',"Multi-run"='#FFBB22')) +
  scale_linetype_manual(values=c("Pivotal" = 'solid', "Subsketching" = 'longdash', "Bootstrap" = 'twodash', "Plug-in" = 'dashed',"Multi-run"='dotted'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_bw()+
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(color = "grey20", size = 15)
  ) +
  xlab('m') +
  ylab('coverage')


ggplot(df_len, aes(x = x)) +
  #ylim(0, 1) +
  geom_point(aes(y = y1, color = "Pivotal",shape='Pivotal'), size = 2) +
  geom_line(aes(y = y1, color = "Pivotal",linetype = 'Pivotal')) +
  geom_point(aes(y = y2, color = "Subsketching",shape="Subsketching"), size = 2) +
  geom_line(aes(y = y2, color = "Subsketching",linetype = 'Subsketching')) +
  geom_point(aes(y = y3, color = "Bootstrap",shape="Bootstrap"), size = 2) +
  geom_line(aes(y = y3, color = "Bootstrap",linetype = 'Bootstrap')) +
  geom_point(aes(y = y4, color = "Plug-in",shape='Plug-in'), size = 2) +
  geom_line(aes(y = y4, color = "Plug-in",linetype='Plug-in')) +
  geom_point(aes(y = y5, color = "Multi-run",shape='Multi-run'), size = 2) +
  geom_line(aes(y = y5, color = "Multi-run", linetype = 'Multi-run')) +
  geom_errorbar(aes(ymin = low[1,], ymax = upp[1,],color = 'Pivotal'), width = 0.05*x[1]) +
  geom_errorbar(aes(ymin = low[2,], ymax = upp[2,],color='Subsketching'), width = 0.05*x[1]) +
  geom_errorbar(aes(ymin = low[3,], ymax = upp[3,],color='Bootstrap'), width = 0.05*x[1]) +
  geom_errorbar(aes(ymin = low[4,], ymax = upp[4,],color='Plug-in'), width = 0.05*x[1]) +
  geom_errorbar(aes(ymin = low[5,], ymax = upp[5,],color='Multi-run'), width = 0.05*x[1]) +
  scale_shape_manual(values = c("Pivotal" = 1, "Subsketching" = 2, "Bootstrap" = 3, "Plug-in" = 4,"Multi-run"=5)) +
  scale_color_manual(values = c("Pivotal" = 'black', "Subsketching" = 'red', "Bootstrap" = 'blue', "Plug-in" = '355E3B',"Multi-run"='#FFBB22')) +
  scale_linetype_manual(values=c("Pivotal" = 'solid', "Subsketching" = 'longdash', "Bootstrap" = 'twodash', "Plug-in" = 'dashed',"Multi-run"='dotted'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_bw()+
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(color = "grey20", size = 15)
  ) +
  xlab('m') +
  ylab('length of confidence interval')


