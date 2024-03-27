library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(wilkoxmisc)
library(reshape2)
library(RColorBrewer)
library(Matrix)

table <- read_tsv("emapper_CDS_combined_contigs_all_pan_samples.tidy.txt")
table <- dcast(table, gene_name ~ Sample, value.var = "Count", fill = 0)
table <- as.data.frame(table)
rownames(table) <- table$gene_name
table$gene_name <- NULL
count <- apply(table,1,nnzero)
table$prevalence <- count


table_filter <- table %>%
filter(prevalence >= 74) #Filter reads in less than 10% of samples

table_filter$prevalence <- NULL
table_filter <- t(table_filter)
table_filter <- as.data.frame(table_filter)
table <- dcast(table, Sample ~ gene_name, value.var = "Count", fill = 0)
write_tsv(table,"eggnog_CDS_table_known_gene_binary.cast.txt")

#make the heatmap with sample type as the colorbar
table <- read_tsv("eggnog_CDS_table_known_gene_binary.cast.txt")
rownames(table) <- table$Sample
table$Sample <- NULL
matrix <- as.matrix(table)

#prepare the color for colorbar and main heatmap
meta <- read_tsv("pan_metadata.txt") %>% filter(!Sample %in% c("SL342522","SL342523")) #remove two samples without any gene identified
meta <- meta %>% arrange(Sample) #reorder the samples
#must be TRUE for all comparsions
my_group <- as.factor(meta$Type)
meta$Type <- factor(meta$Type, levels = c("Skin","Indoor surface","Subway air","Subway surface", "Urban public surface","Pier surface"))
colSide <- brewer.pal(8,"Paired")[my_group]
#colMain <- colorRampPalette(brewer.pal(8, "Set1"))(2)
colMain <- c("#e0f3f8","#4575b4")

heatmap <- heatmap(matrix, scale="none", xlab="Gene", ylab="Sample", cexRow=0.1, col=colMain,RowSideColors=colSide,labRow=NULL, labCol=NULL)
quartz.save(type="png",dpi=1200,file="eggnog_heatmap_by_contig_binary_w_bin_name.png")

heatmap <- heatmap(matrix, scale="none", xlab="Gene", ylab="Sample", margins = c(2, 2), col=colMain,RowSideColors=colSide,labRow=NA, labCol=NA)
quartz.save(type="png",dpi=1200,file="eggnog_heatmap_by_contig_binary_final.png")


