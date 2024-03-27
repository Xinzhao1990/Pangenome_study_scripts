library(readr)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
    
table <- read_tsv("mluteus_accessory_core_BGCs.txt")

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')), space='Lab')

library(tidyverse)
library(grid)

hacky_df <- data.frame(
var = c("Core", "Global BE", "HK BE", "NCBI"),
var_color = c("#e41a1c", "#377eb8","#4daf4a", "#984ea3")
)
# plot code
plot_new <-
ggplot(df) +    # don't specify x and y here.  Otherwise geom_rect will complain.
geom_rect(
data=hacky_df,
aes(xmin=-Inf, xmax=Inf, fill=var_color, alpha=0.4)) +
geom_point(aes(x = Genome, y = Gene_cluster_type)) +
coord_cartesian(clip="off") +
facet_wrap(~ Source, scales = "free") +
scale_fill_manual(values = c("Core" = "#e41a1c", "Global BE" = "#377eb8", "HK BE" = "#4daf4a", "NCBI" = "#984ea3")) +
theme_bw() +
theme(
strip.background = element_rect(fill=NA),
strip.text = element_text(face="bold")
)

plot_new

table2 <- read_tsv("mluteus_accessory_core_BGCs.txt")
Plot <- ggplot(table2,aes(x=Genome, y=Gene_cluster_type))
Plot <- Plot + facet_grid(~Source, scales = "free",space="free") #+ scale_y_continuous(expand = c(0,0))
#Plot <- Plot + geom_point(aes(size=read_coverage,color=Sample_type), alpha=0.5)
Plot <- Plot + geom_point(aes(fill=Gene_cluster_type,size=3),pch=21,alpha = 0.8)
Plot <- Plot + scale_fill_manual(values = c("Core" = "#e41a1c", "Global BE" = "#377eb8", "HK BE" = "#4daf4a", "NCBI" = "#984ea3"))
Plot <- Plot + theme(strip.background = element_rect(fill = "Dark2"))
#Plot <- Plot + theme_bw()

Plot <- Plot + theme(axis.text.x = element_blank())
Plot <- Plot + theme(axis.ticks.x = element_blank())
#Plot <- Plot + theme(axis.text.x = element_text(angle = 90, hjust = 0))
Plot <- Plot + scale_fill_brewer(palette = "Set1")
Plot <- Plot + xlab(paste0("M.luteus genomes")) + ylab(paste0("Type of BGCs"))
Plot <- Plot + theme(legend.position = "bottom")



table <- read_tsv("identical_MAGs_summary_meta.txt") %>%
    select(Genotype,Taxonomy,Type, Location)
table$count <- 1

table <- table %>%
    group_by(Genotype,Taxonomy,Type, Location) %>%
    mutate(sum=sum(count)) %>%
    ungroup() %>%
    select(-count) %>%
    unique()
# identify idential genome across occupants and residence
Plot <- ggplot(table,aes(x=Location, y=Taxonomy, shape=Genotype))
Plot <- Plot + facet_wrap(~Type)
#Plot <- Plot + geom_point(aes(size=read_coverage,color=Sample_type), alpha=0.5)
Plot <- Plot + geom_point(aes(color=Location,size=sum),alpha = 0.7)
#Plot <- Plot + scale_shape_manual(values = c(16,17,15,3,7,8,9))
Plot <- Plot + theme_bw()
Plot <- Plot + theme(axis.text.x = element_text(angle = 90, hjust = 0))
Plot <- Plot + scale_color_brewer(palette = "Set2")
#Plot <- Plot + theme(axis.title = element_text(size=16))
#Plot <- Plot + theme(legend.title = element_text(size=10))
#Plot <- Plot + theme(legend.text = element_text(size=8))
Plot <- Plot + ylab(paste0("Taxonomy of representative MAG"))+xlab(paste0("Residence"))
