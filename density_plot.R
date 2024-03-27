library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

OTU <- read_tsv("kraken_species_abundance_w_classified.tidy.txt") %>%
    select(-Count)
Meta <- read_tsv("pan_metadata.txt") %>% select(Sample, Type)

Merge <- left_join(OTU,Meta)

#dominant species (mean RA > 1%)
Merge <- Merge %>% filter(Species %in% c("Micrococcus luteus","Cutibacterium acnes","Gordonia bronchialis","Moraxella osloensis", "Janibacter indicus", "Dermacoccus nishinomiyaensis", "Kocuria palustris", "Xanthomonas campestris","Bradyrhizobium sp. BTAi1","Kytococcus sedentarius"))

#indicator species
Merge <- Merge %>% filter(Species %in% c("Vibrio alginolyticus","Gloeocapsa sp. PCC 7428","Staphylococcus capitis"))

Merge$Type <- factor(Merge$Type, levels = c("Skin","Indoor surface","Subway air","Subway surface", "Urban public surface","Pier surface"))
Plot <- ggplot(Merge,aes(x=RelativeAbundance,colour=Type)) + geom_density(lwd = 1) + facet_wrap(~Species,scales="free")
Plot <- Plot + theme_bw()
Plot <- Plot + xlab(paste0("Relative Abundance (%)")) + ylab(paste0("Density")) + theme(legend.position="bottom")
Plot <- Plot + scale_color_brewer(palette = "Paired")