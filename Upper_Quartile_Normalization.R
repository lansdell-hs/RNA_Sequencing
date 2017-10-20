setwd("C:/Users/hslansd/Desktop")

library(factoextra)
library(FactoMineR)
library(caret)
library(plyr)

1RSEMreads<-read.csv("rsem_raw_batch1.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
2RSEMreads<-read.csv("rsem_raw_batch2.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

data<-as.matrix(rbind(t(1RSEMreads),t(2RSEMreads)))

rownames(data) <- substring(rownames(data), 2)

data<-data[,colSums(data)>0]

data<-t(data)
data<-as.data.frame(data)

seventyfive_quart<-numcolwise(function(x) quantile(x, .75))(data)


for(i in 1:192){data[,i]<-data[,i]/seventyfive_quart[,i]}  #192 number of columns


mean_quart<-sum(seventyfive_quart)/192

for(i in 1:192){data[,i]<-data[,i]*mean_quart}

data<-t(data)


Ext_list=c("XXXXXX","XXXX","XXXX")

data<-data[!rownames(data) %in% Ext_list,]
fviz_pca_ind(prcomp(data))
