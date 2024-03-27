library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(wilkoxmisc)
library(reshape2)
library(RColorBrewer)

table <- read_tsv("selected_kegg_ko_module_biodagradation_anti_resistance.cast.txt")
table$Module_step <- NULL
table$KO <- NULL
table$Module <- NULL
table$Module_Category <- NULL

meta <- read_tsv("selected_kegg_ko_module_biodagradation_anti_resistance.cast.txt") %>%
    select(KO2,Module_Category)

list <- read_tsv("unique_ko_module_genome_list.txt")
tax <- read_tsv("pan_gtdbtk_summary_meta_checkm_860.tsv") %>%
    select(Genome,Phylum)

list <- left_join(list,tax)
#make the heatmap with sample type as the colorbar
rownames(table) <- table$KO2
table$KO2 <- NULL
rownames(table) == meta$KO2 #must be TRUE for all
colnames(table) == list$Genome #must be TRUE for all
func <- as.factor(meta$Module_Category)
rowSide <- brewer.pal(8,"Dark2")[func]
phylum_col <- as.factor(list$Phylum)
colorset <- colorRampPalette(brewer.pal(12,"Paired"))(12)
colSide <- colorset[phylum_col]
matrix <- as.matrix(table)

my_group <- as.factor(colnames(table))

#colMain <- colorRampPalette(brewer.pal(8, "Purples"))(25)
colMain <- c("black","grey")

library(heatmap3)
heatmap <- heatmap3(matrix, scale="none", Rowv=NA, labCol=NA,cexRow=0.3,margins = c(1, 5),col=colMain,RowSideColors=rowSide)

legend(x="top", legend=c("Aromatics degradation", "Drug efflux transporter/pump", "Drug resistance"),fill=colorRampPalette(brewer.pal(12,"Paired"))(13))

heatmap <- heatmap3(matrix, scale="none", Rowv=NA, Colv=NA, labCol=NA,cexRow=0.3,margins = c(1, 5),col=colMain,ColSideColors=colSide,RowSideColors=rowSide)

legend(x="top", legend=c("Acidobacteriota","Actinobacteriota", "Bacteroidota", "Chloroflexota", "Cyanobacteria", "Deinococcota",     "Eremiobacterota", "Firmicutes", "Gemmatimonadota", "Myxococcota", "Planctomycetota", "Proteobacteria"), fill=colorRampPalette(brewer.pal(12,"Paired"))(12))


quartz.save(type="pdf",dpi=1200,file="unique_selected_ko_module_cluster_A_NEW.pdf")


