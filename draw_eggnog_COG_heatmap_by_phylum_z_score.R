library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(wilkoxmisc)
library(reshape2)
library(RColorBrewer)

table <- read_tsv("combined_COG_all_genomes_COG_percentage.txt")
genome <- read_tsv("m_luteus_genome_metadata.txt")

names(table)[1] <- c("Genome")

merge <- left_join(table,genome) %>%
    filter(!COG_category == "S")

#calculate the z-score for making the heatmap

table <- dcast(merge, Genome2 ~ COG_category, value.var = "Percentage", fill = 0)

#make the heatmap with sample type as the colorbar
genome <- genome %>%
arrange(Genome2)

table$Genome2 == genome$Genome2
rownames(table) <- genome$Genome2
table$Genome2 <- NULL
matrix <- as.matrix(table)

#prepare the color for colorbar and main heatmap
#meta <- read_tsv("pier_metadata_complete.txt") %>%
    filter(Sample %in% rownames(table))
rownames(matrix) == genome$Genome2 #must be TRUE for all comparsions
my_group <- as.factor(genome$Source)
colSide <- brewer.pal(8,"Set1")[my_group]
colMain <- colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(50)
#colMain <- colorRampPalette(brewer.pal(8, "Purples"))(50)

#make heatmap
heatmap <- heatmap(matrix, scale="none", xlab="COG category", margins = c(3, 6), col=colMain, RowSideColors=colSide)

legend(x="topleft", legend=c("Core pangenome", "Global BE","HK BE", "NCBI"),fill=colorRampPalette(brewer.pal(8, "Set1"))(8))

quartz.save(type="pdf",dpi=1200,file="eggnog_COG_percentage_mluteus_heatmap_RdYlBu.pdf")

library(heatmap3)
heatmap <- heatmap3(matrix, scale="none", ylab="Sample-based contig", xlab="COG category",margins = c(3, 6), col=colMain,RowSideColors=colSide)
quartz.save(type="pdf",dpi=1200,file="eggnog_COG_percentage_mluteus_heatmap_legend_RdYlBu.pdf")

