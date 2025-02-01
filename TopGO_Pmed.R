#Install Biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Install topGO
BiocManager::install(c("topGO"))
library(topGO)
setwd("")
#loading your data
#Reading in GO annotations for the genes
gene2GO_map <- readMappings(file="Gene2GO_pmed.txt", sep = "\t")
#Defining your list of genes of interest, and the 'gene universe' 
geneUniverse<-names(gene2GO_map)
head(geneUniverse)
genesOfInterest<- read.table("72_down genes.txt",header=FALSE) #Genes you are interested in 
genesOfInterest <- as.character(genesOfInterest$V1)

#Tell TopGO where these interesting genes appear in the 'geneUniverse' vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#Putting the data together into an R object
myGOdata <- new("topGOdata",
                description="My project",
                ontology= "BP", 
                allGenes=geneList, 
                annot = annFUN.gene2GO, gene2GO = gene2GO_map)
myGOdata

# The list of genes of interest can be accessed using the method sigGenes():
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)

# Performing enrichment tests (Fisher's exact test based on gene counts)
resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
#p.vals_res <- score(resultFisher)
#adj.p.vals_res <- p.adjust(p.vals_res)
resultFisher

#The p-values have not been corrected for multiple testing.
#We can list the top ten significant results found:
allRes <- GenTable(myGOdata, 
                   classicFisher = resultFisher, 
                   orderBy = "resultFisher", 
                   ranksOf = "classicFisher", 
                   topNodes = 50)
allRes
#p.vals <- as.vector(allRes[,6]); # ONLY DO WHATS IS BELOW IF YOU WANT TO CORRECT FOR MULTIPLE TESTING (ADJUST P-VALUE). 
#NOT NECCESSARY PLS SEE:https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/

#classicFisher.adj_p <- as.vector(adj.p.vals_res[allRes[,1]]);
#newRes <- cbind(allRes, classicFisher.adj_p);

#newRes

#sel.terms <- as.vector(newRes[,1]);
#sel.terms;
#num.ann.genes <- countGenesInTerm(GOdata, sel.terms);
#ann.genes <- genesInTerm(GOdata, sel.terms);
#ans = c("GO","geneID");

#for (i in 1:length(ann.genes)) {for (j in ann.genes[i]) {for (k in j) {if (k %in% as.vector(u_myInterestedGenes)) ans=rbind(ans,c(names(ann.genes[i]),(k))) }}};

write.table(allRes, file="72_down_genes.csv", sep=",")
#We can visualise the position of the statistically significant GO terms in the GO hierarchy by using the following functions:
install.packages("BiocManager")
BiocManager::install("Rgraphviz")
library("Rgraphviz")
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "72down_classi_5_all", useInfo = "all", pdfSW = TRUE)

#Finding the genes annotated with significant GO terms
length(usedGO(myGOdata))
myterms = c("GO:0006928", "GO:0007018", "GO:0015979")
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}

# above method came from http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html