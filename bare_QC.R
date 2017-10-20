setwd("C:/Users/hslansd/Desktop/Total RNA/")
totalRNA<-read.csv("data.csv", header = TRUE, row.names =1, stringsAsFactors = FALSE )

data<-as.matrix(t(totalRNA))
rownames(data) <- substring(rownames(data), 2)
# cols at this point should be genes n= 20530, rows are samples n= 176
data<-data[,-nearZeroVar(data)]
data<-data<-data[,colSums(data)>0]
#cols genes n=20353 row samples n=176
Ext_list=c("XXXX","XXXX","XXX") #removed for privacy

fviz_pca_ind(prcomp(data))

data<-t(data)
plot(hclust(as.dist(1-cor(data))))
data<-data[,!colnames(data) %in% Ext_list]

var_genes <- apply(data, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]

highly_variable <- data[select_var,]
dim(highly_variable)
head(highly_variable)
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlOrBr")
morecols <- colorRampPalette(mypalette)
library(gplots)

#heatmap.2(highly_variable,col=rev(morecols(50)),trace="none", main="Top 100 most variable genes w/outliers marked (cluster)", scale = "row",ColSideColors=ifelse(colnames(highly_variable) %in% Ext_list,"red","blue"))
heatmap.2(highly_variable,col=rev(morecols(50)),trace="none", main="Top 100 most variable genes across samples",scale="row")

write.csv(file ="TotalRNA_postQC.csv" ,data)
