library(readstata13)
library(dplyr)
library(caret)

usa1940 <- readstata13::read.dta13("D:/Dropbox/R/sketching/code_randomized_algorithms/large real data/usa1940.dta")
df0<-filter(usa1940,ind=='0124') 
#Industry code Educational service "https://usa.ipums.org/usa/volii/ind1940.shtml"

df1<-droplevels(df0)
df2<-subset(df1,select=-c(ind,ln_hr_wage))

dummy <- dummyVars(" ~ .", data=df2,fullRank = TRUE)

df3 <- data.frame(predict(dummy, newdata=df2))
df4<-df3[,-(which(colSums(df3)==1))]

occ_count <- sum(grepl("occ", names(df4)))
n<-nrow(df4);p<-ncol(df4)
c<-rep(0,p)
c[p-occ_count]<-1

y<-as.vector(df0$ln_hr_wage)
X<-as.matrix(df4)
