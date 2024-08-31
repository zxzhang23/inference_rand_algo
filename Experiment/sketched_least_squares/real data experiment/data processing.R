library(readstata13)
library(dplyr)
library(caret)

usa1940 <- readstata13::read.dta13("usa1940.dta")
summary(usa1940)
df0<-filter(usa1940,ind=='0124') #Industry code Educational service "https://usa.ipums.org/usa/volii/ind1940.shtml"
df1<-subset(df0,select=-c(bpl,occ,ind))

dummy <- dummyVars(" ~ .", data=df1,fullRank = TRUE)

#perform one-hot encoding on data frame
df3 <- data.frame(predict(dummy, newdata=df1))
df<-df3[,-(which(colSums(df3)==0))]

save(df,file="usa1940edu.rda")
