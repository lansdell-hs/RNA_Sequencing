source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")

library(DESeq) 
library(statmod) 
library(pcaMethods)
library(fastICA)
require(DESeq)

1mRNA<- read.csv("batch1.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
2mRNA<- read.csv("batch2.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

1mRNA_mat<-as.matrix(f96mRNA)
2mRNA_mat<-as.matrix(s96mRNA)

data <- rbind(1mRNA_mat, 2mRNA_mat)
data<-data<-data[rowSums(data)>0,]


means <- rowMeans(data)
vars <- apply(data,1,var)
cv2 <- vars/means^2
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))

require(statmod)
minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= minMeanForFit 
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(data) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- data[varorder,]

par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 1000 genes
points(log(means[varorder[1:1000]]),log(cv2[varorder[1:1000]]),col=2)

m <- oed[1:1000,]
heatmap(m/apply(m,1,max),zlim=c(0,1),col=topo.colors(50),Rowv=NA,Colv=NA,labRow=NA,scale="column")

pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
table(sigVariedGenes)
