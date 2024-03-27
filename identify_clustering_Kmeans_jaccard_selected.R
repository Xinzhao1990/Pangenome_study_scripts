#Generate bray-curtis and jaccard
library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)
library(vegan)
library(reshape2)
library(RVAideMemoire)
library(ggplot2)
library(fpc)
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

tax_cast <- read_tsv("combined_all_genomes_kegg_ko_unique_binary.cast.txt")
names(tax_cast)[1] <- c("Genome")

meta <- read_tsv("pan_gtdbtk_summary_meta_checkm_1663.tsv")
meta$Genome <- gsub("_bin","bin",meta$Genome)
meta <- meta %>%
filter(Completeness >=90)

tax_cast2 <- tax_cast %>% filter(Genome %in% meta$Genome)
write_tsv(tax_cast2,"combined_hq_genomes_kegg_ko_unique_binary.cast.txt")

tax_cast <- read.table("combined_all_genomes_kegg_ko_unique_binary_selected.cast.txt",header=T,row.names=1)
#tax_cast <- t(tax_cast)
#distance_b <- vegdist(tax_cast, method="bray")
distance_j <- vegdist(tax_cast, method="jaccard")

BrayMatrix <- as.matrix(distance_j)

library(factoextra)
library(cluster)
#Dertemines and visualize the optimal number of clusters
fviz_nbclust(BrayMatrix, kmeans, method="wss") #for kmeans
fviz_nbclust(BrayMatrix, pam, method="silhouette") #for pam

#Prediction strength
ps <- prediction.strength(BrayMatrix,M=100)
# print the prediction strength values
#Prediction strength
#Clustering method:  kmeans
#Maximum number of clusters:  10
#Resampled data sets:  100
#Mean pred.str. for numbers of clusters:  1 0.969107 0.8396755 0.7482511 0.53613 0.4916278 0.4602628 0.4354411 0.4195738 0.4135019
#Cutoff value:  0.8
#Largest number of clusters better than cutoff:  3

Prediction strength
Clustering method:  kmeans
Maximum number of clusters:  10
Resampled data sets:  100
Mean pred.str. for numbers of clusters:  1 0.9101298 0.7041357 0.6813059 0.5656991 0.5034666 0.4842946 0.4459961 0.4130293 0.3983564
Cutoff value:  0.8
Largest number of clusters better than cutoff:  2

#Calculate number of clusters and their clustering strength (cluster package)
library(cluster)
#cluster <- pam(BrayMatrix,2,cluster.only=TRUE)
pm.res <- pam(BrayMatrix,6)
plot(pm.res)
km.res <- kmeans(BrayMatrix, 2, nstart = 50) #random starts do help here with too many clusters eg cluster=5
cluster <- km.res$cluster
write.table(cluster,"kmeans_2_cluter_binary_test2024.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#Calculate number of clusters using k-mean clustering (stats package)
km <- kmeans(BrayMatrix, 2, nstart = 50)
cluster <- km$cluster
write.table(cluster,"kmeans_cluster_2.txt", sep="\t", col.names=TRUE, row.names=TRUE)


#make principle component plot
my_data_withpam = cbind(BrayMatrix, cluster=km.res$cluster,shape=meta$Type)
plot <- fviz_cluster(km.res,data=BrayMatrix,palette = "Accent", ellipse.type="t", #Concentration ellipse
repel = TRUE,#Avoid label overplotting #data=BrayMatrix is required for k means
show.clust.cent=FALSE,
geom = "point",
ggtheme = theme_classic())
#plot <- plot + geom_point(aes(shape=meta$Type))
plot <- plot + ggtitle("PAM clustering") + theme(plot.title = element_text(hjust = 0.5))
#fviz_cluster plot is equivalent to this one
source("ggbiplot.r")
table.pca <- prcomp(BrayMatrix,center = TRUE,scale. = TRUE)
km.res$cluster <- gsub("1","A",km.res$cluster)
km.res$cluster <- gsub("2","B",km.res$cluster)
plot <- ggbiplot(table.pca,group=km.res$cluster,ellipse=TRUE,var.axes=FALSE)
plot <- plot + scale_color_manual(name="Cluster", values=c("orange", "purple"))
plot <- plot + scale_shape_manual(name="Type", values=c(1:3,17:19))
plot <- plot + geom_point(size=3, aes(color=km.res$cluster,shape=meta$Type))
plot <- plot + scale_color_brewer(palette = "Set1")

#make PCA plot using ggplot2
source("ggbiplot.r")
table.pca <- prcomp(BrayMatrix,center = TRUE,scale. = TRUE)
km.res$cluster <- gsub("1","A",km.res$cluster)
km.res$cluster <- gsub("2","B",km.res$cluster)
km.res$cluster <- gsub("3","C",km.res$cluster)
km.res$cluster <- gsub("4","D",km.res$cluster)
km.res$cluster <- gsub("5","E",km.res$cluster)
km.res$cluster <- gsub("6","F",km.res$cluster)
table.pca.df <- data.frame(summary(table.pca)$x)
xlab <- paste0("PC1 (Variance Explained: 36.1%)")
ylab <- paste0("PC2 (Variance Explained: 24.5%)")

meta <- read_tsv("pan_gtdbtk_summary_meta_checkm_1663.tsv")
meta$Genome <- gsub("_","",meta$Genome)
meta <- meta %>%
filter(meta$Genome %in% rownames(BrayMatrix))
meta <- meta %>% arrange(meta$Genome)
rownames(BrayMatrix) == meta$Genome #MUST be TRUE

library(RColorBrewer)

#colorCount = length(unique(meta$Phylum))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))

meta$Phylum <- gsub("p__","", meta$Phylum)
meta$Class <- gsub("c__"," ", meta$Class)

km.res$cluster <- factor(km.res$cluster, levels = c("A", "B"))

Plot <- ggplot(table.pca.df, aes(x = PC1, y = PC2, shape=km.res$cluster))
Plot <- Plot + geom_point(size=3.5,alpha=0.8,aes(fill=meta$Phylum)) + stat_ellipse(type = "norm", linewidth = 1.5, linetype = 2, aes(color = km.res$cluster))+stat_ellipse(type = "t",aes(color = km.res$cluster)) #add pch=21 to show legend
#Plot <- stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t")
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_bw()
Plot <- Plot + scale_shape_manual(values=c(21,23))#, 24,25))#, 24,25,4))
Plot <- Plot + scale_fill_brewer(palette = "Set1")
Plot <- Plot + scale_color_brewer(palette = "Dark2")
#Plot <- Plot + scale_color_manual(values= getPalette(colorCount))
Plot <- Plot + theme(axis.title=element_text(size=14))
Plot <- Plot + theme(axis.text=element_text(size=10))
Plot <- Plot + ggtitle("K-means clustering (k=2)")+theme(plot.title = element_text(hjust = 0.5))


colorCount = length(unique(meta$Class))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

Plot <- ggplot(table.pca.df, aes(x = PC1, y = PC2, shape=km.res$cluster))
Plot <- Plot + geom_point(size=2,alpha=0.8,aes(fill=meta$Class)) + stat_ellipse(type = "norm", linetype = 2, aes(color = km.res$cluster))+stat_ellipse(type = "t",aes(color = km.res$cluster)) #add pch=21 to show legend
#Plot <- stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t")
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_bw()
Plot <- Plot + scale_shape_manual(values=c(21,23,24))#, 24,25))#, 24,25,4))
#Plot <- Plot + scale_fill_brewer(palette = "Set1")
Plot <- Plot + scale_color_brewer(palette = "Accent")
Plot <- Plot + scale_fill_manual(values= getPalette(colorCount))
Plot <- Plot + theme(axis.title=element_text(size=12))
Plot <- Plot + theme(axis.text=element_text(size=10))
Plot <- Plot + ggtitle("K-means clustering (k=3)")+theme(plot.title = element_text(hjust = 0.5))


