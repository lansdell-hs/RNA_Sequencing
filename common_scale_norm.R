setwd("C:/Users/hslansd/Desktop")
library(factoextra)
library(FactoMineR)
library(caret)

1RSEMreads<-read.csv("rsem_rawbatch1.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
2RSEMreads<-read.csv("rsem_rawbatch2.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

data<-as.matrix(rbind(t(1RSEMreads),t(2RSEMreads)))

rownames(data) <- substring(rownames(data), 2)


fviz_pca_ind(prcomp(data))

data<-data[,-nearZeroVar(data)]
data<-data<-data[,colSums(data)>0]


data<-t(data)
plot(hclust(as.dist(1-cor(data))))

dist_val<-1-cor(data)


library(pracma)
pseudoRefSample<-geomean(data)
nz<-pseudoRefSample>0
ratios<- data[,nz]/pseudoRefSample[nz]

sizeFactors<-0
for(i in 1:1751){sizeFactors[i]<-median(ratios[,i])}

normCounts<-data[,nz]/sizeFactors


normCounts<-t(normCounts)
plot(hclust(as.dist(1-cor(normCounts))))
plot(normCounts)
