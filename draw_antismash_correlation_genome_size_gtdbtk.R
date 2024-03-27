library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

table <- read_tsv("antismash_BGC_count_per_bin.txt") %>%
    select(-Gene_cluster_type) %>%
    group_by(Bin) %>%
    mutate(total =sum(number)) %>%
    select(-number) %>%
    unique()
size <- read_tsv("genome_size_860_MAGs.txt")
merge <- left_join(table,size)
genome <- read_tsv("pan_gtdbtk_summary_meta_checkm_1663.tsv")
names(genome)[1] <- c("Bin")
#genome$Bin <- gsub("bin.","bin", genome$Bin)
merge <- left_join(merge,genome)

merge$Phylum <- gsub("p__","",merge$Phylum)
#list <- read.csv("Wdb_100_dereplicated_genomes.csv")
#merge <- merge %>% filter(Bin %in% list$genome)

#calculate the coefficient for regression line passing the origin
lm <- lm(total ~ 0 +genome_size_MB,data=merge)

Call:
lm(formula = total ~ 0 + genome_size_MB, data = merge)

Coefficients:
genome_size_MB
1.523

formula <- merge$total ~ 0 + merge$genome_size_MB

merge2 <- read_tsv("antismash_BGC_count_genome_size_taxonomy.txt")
test <- cor.test(merge2$total,merge2$genome_size_MB,method="spearman")

Spearman's rank correlation rho

data:  merge2$total and merge2$genome_size_MB
S = 51421831, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho
0.5149304

> test$p.value
[1] 2.090848e-59

library(RColorBrewer)

colorCount = length(unique(merge2$Phylum))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

merge2 <- merge2 %>% select(genome_size_MB,total,Phylum)
merge2 <- as.data.frame(merge2)
plot <- ggplot(merge2,aes(x=genome_size_MB,y=total,color=Phylum))
plot <- plot + geom_point(size=2, alpha=0.8) #+ geom_smooth(method='loess',color="black")
plot <- plot + theme_test()
plot <- plot + scale_color_manual(values= getPalette(colorCount))
#plot <- plot + scale_x_discrete(breaks = seq(0, 9, by = 1))
#plot <- plot + abline(a=0,b=1.68)
plot <- plot + xlab(paste0("Genome size (Mb)")) + ylab(paste0("Number of BGC"))

merge2 <- merge %>%
    filter(!Phylum %in% c("p__Firmicutes","p__Gemmatimonadetes"))

library(wilkoxmisc)
kruskal.test(total~Phylum, data=merge2)

Kruskal-Wallis rank sum test

data:  total by Phylum
Kruskal-Wallis chi-squared = 44.833, df = 5, p-value = 1.569e-08

library(pgirmess)

kruskalmc(total~Phylum, data=merge2)

merge$Bin <- gsub("_"," ",merge$Bin)
merge$Sample <- word(merge$Bin,1)

meta <- read_tsv("pier_metadata_complete.txt") %>%
	select(Sample,Surface_material,Sample_type,Group, Location)
	
merge <- left_join(merge,meta)

f__Chroococcidiopsidaceae <- merge %>%
    filter(Family == "f__Chroococcidiopsidaceae") %>%
    filter(!Gene_cluster_type == "other")
Plot <- ggplot(f__Chroococcidiopsidaceae,aes(x=Bin, y=number,fill=Gene_cluster_type))
Plot <- Plot + geom_bar(stat="identity",width=0.5)
#Plot <- Plot + facet_grid(Sample_type~Location,scales = "free_y")+ scale_y_continuous(expand = c(0,0))
Plot <- Plot + facet_wrap(~Species,scales = "free_y",nrow=1)+ scale_y_continuous(expand = c(0,0))
Plot <- Plot + theme(axis.text.x = element_text(angle=-90, size=8))
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + xlab(paste0("Erythrobacteraceae genome")) + ylab(paste0("BGC count"))
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme(axis.title = element_text(size=16))
Plot <- Plot + theme(legend.title = element_text(size=12))
Plot <- Plot + theme(legend.text = element_text(size=10))



s__Deinococcus_swuensis <- merge %>% filter(Species == "s__Deinococcus_swuensis")
Plot <- ggplot(s__Deinococcus_swuensis ,aes(x=Bin, y=number,fill=Gene_cluster_type))
Plot <- Plot + geom_bar(stat="identity",width=0.5)
Plot <- Plot + facet_grid(Sample_type~Location,scales = "free_y")+ scale_y_continuous(expand = c(0,0))
#Plot <- Plot + facet_wrap(~Location,scales = "free_y",nrow=1)+ scale_y_continuous(expand = c(0,0))
Plot <- Plot + theme(axis.text.x = element_blank())
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + xlab(paste0("Deinococcus swuensis genome")) + ylab(paste0("BGC count"))
Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title = element_text(size=16))
Plot <- Plot + theme(legend.title = element_text(size=12))
Plot <- Plot + theme(legend.text = element_text(size=10))
