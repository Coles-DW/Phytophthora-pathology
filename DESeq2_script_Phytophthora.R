####DESEQ2 OZMYC data############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
##################################
library("DESeq2")
library("ggplot2")
library('reshape')
workingDir = ""
setwd(workingDir)
data=read.table("filtered_unique_counts.txt", header = TRUE, row.names = 1)
colnames(data)
sample=colnames(data)
condition = c("mock","mock","mock","mock","time12","time12","time12","time12","time24","time24",
              "time24","time24","time72","time72","time72","time72")
colData=cbind(colnames(data), sample, condition)
######verification des col -treatment###
colData=data.frame(colData)
write.table(colData, file="PmedCol_Samples_Conditions.txt", sep="\t")
######boxplot of data transnformation
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(data + epsilon)), breaks=100, col="blue", border="white",
     main="Expressed_Log2-transformed counts per gene_Pmed", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
dev.copy(png,'Log2-transformed_Histogram_Pmed.png')
dev.off()
condition_1 <- as.data.frame(t(condition))
boxplot(log2(data + epsilon), condition_1$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)")
dev.copy(png,'boxplot_expressed_transformeddata_Pmed.png')
dev.off()
######dds=data deseq
dds= DESeqDataSetFromMatrix(data, colData, ~condition)
dds= DESeq(dds)
results(dds)
str(dds)
######sortie donnees normalises############
norm_counts=counts(dds, normalized=TRUE)
write.table(norm_counts, file="Pmed_normalized_expressed.txt", sep="\t")
##rlog tranformation
rld<- rlogTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)

library("gplots")
library('RColorBrewer')
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(10, 10))
dev.copy(png,'Expressed_Pmed_sampletosample_heatmap.png')
dev.off()
print(plotPCA(rld))
pcaData <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=1) + 
  geom_mark_ellipse(aes(fill = condition,
                        colour = condition)) +
  theme(legend.position = 'bottom') +
  coord_equal() +
  xlim(-100,100)+
  ylim(-50,50) +
  xlab(paste0("\n PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance \n")) + 
  theme_classic()
p

dev.copy(pdf, 'Expressed_Pmed_PCA2.pdf')
dev.off()
#####savoir la dimension du fichier, nombre col et lignes)###
dim(norm_counts)
str(dds)
#####Filtering 10 reads per sample repetition->excel########
##############Hclust##########
## transfo en log2(x+1) ##
norm_counts_log=log2(norm_counts+1)
## calcul correlation de Pearson##
r = cor(norm_counts_log, method = "pearson")
write.table(r, file="Pmed_correlation_values.txt", sep="\t")
## calcul de la dissimilarité ##
d=1-r
## Clustering ##
h = hclust(as.dist(d), method = "complete", members = NULL)
plot(h)
dev.copy(pdf, 'Expressed_Pmed_hclust.pdf')
dev.off()

####comparaison####
####Pairwise comparison#####
res1 <- results(dds, contrast=c("condition","time12","mock"))
res2 <- results(dds, contrast=c("condition","time24","mock"))
res3 <- results(dds, contrast=c("condition","time72","mock"))


######Mise en forme#########
res1=res1[,c(2,6)]
res2=res2[,c(2,6)]
res3=res3[,c(2,6)]


#Fichier Final
colnames(res1)=c("Log2Foldchange_12hourvsmock","padj_12hourvsmock")
colnames(res2)=c("Log2Foldchange_24hourvsmock","padj_24hourvsmock")
colnames(res3)=c("Log2Foldchange_72hourvsmock","padj_72hourvsmock")

chickpeageneexpressionresults=cbind(res1, res2, res3)
write.table(chickpeageneexpressionresults, file="Differentialexpression_Pmed.txt", sep="\t")



