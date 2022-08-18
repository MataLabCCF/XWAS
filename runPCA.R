library(GENESIS)
library(SeqArray)
library(SNPRelate)
library(GWASTools)


options <- commandArgs(trailingOnly = TRUE)

vcfName = options[1]
gdsName = options[2]
csvName = options[3]


snpgdsVCF2GDS(vcfName, gdsName,  method="biallelic.only")
genofile <- openfn.gds(gdsName)
pca<-snpgdsPCA(genofile,  autosome.only=FALSE)

tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],    
                  PC4 = pca$eigenvect[,4],    
                  PC5 = pca$eigenvect[,5],
                  PC6 = pca$eigenvect[,6],
                  PC7 = pca$eigenvect[,7],
                  PC8 = pca$eigenvect[,8],
                  PC9 = pca$eigenvect[,9],
                  PC10 = pca$eigenvect[,10],
                  stringsAsFactors = FALSE)

write.table(tab, csvName, row.names=FALSE, quote=FALSE, sep = "\t") 
