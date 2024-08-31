library(readstata13)
library(dplyr)
library(caret)

usa1940 <- readstata13::read.dta13("usa1940.dta")
summary(usa1940)
df0<-filter(usa1940,ind=='0124') #Industry code Educational service "https://usa.ipums.org/usa/volii/ind1940.shtml"
df1<-subset(df0,select=-c(bpl,occ,ind))

#perform one-hot encoding on dataframe
dummy <- dummyVars(" ~ .", data=df1,fullRank = TRUE)
df3 <- data.frame(predict(dummy, newdata=df1))
df<-df3[,-(which(colSums(df3)==0))]

save(df,file="usa1940edu.rda")


######################################################################
load("D:/Dropbox/R/sketching/code_randomized_algorithms/large real data/usa1940edu.rda")

y<-as.numeric(df$ln_hr_wage)
X<-as.matrix(subset(df,select=-c(ln_hr_wage)))
