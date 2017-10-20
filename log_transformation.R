setwd("C:/Users/hslansd/Desktop")
library(factoextra)
library(FactoMineR)
library(caret)


1RSEMreads<-read.csv("rsem_rawbatch1.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
2RSEMreads<-read.csv("rsem_rawbatch2.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

data<-as.matrix(rbind(t(1RSEMreads),t(2RSEMreads)))

rownames(data) <- substring(rownames(data), 2)


data<-data[,-nearZeroVar(data)]
data<-data<-data[,colSums(data)>0]

smallNonZeroNum <- 1
data[data==0]<-smallNonZeroNum

data[] <- vapply(data, log2, numeric(1))

data<-t(data)
plot(hclust(as.dist(1-cor(data))))

fviz_pca_ind(prcomp(data))

