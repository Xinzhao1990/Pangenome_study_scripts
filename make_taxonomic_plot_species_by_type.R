library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)

#Open taxonomy OTU table
OTU <- read_tsv("merge_abundance_species_kraken_pan.tidy.txt") %>%
mutate(RelativeAbundance = fraction_total_reads*100)%>%
select(name,Sample,RelativeAbundance)
names(OTU)[1] <- c("Species")
Meta <- read_tsv("pan_metadata.txt") %>% select(Sample,Type)

#Merge relative abundance table and metatable together
OTUTable <- merge(OTU, Meta, by = "Sample", all.x = TRUE)

OTUTable  <- OTUTable %>% filter(Type == "Pier surface")
#collapse taxa table to only 5 or 8 top phyla, genus, family, etc (require reshape2).
OTUTable  <- collapse_taxon_table(subtable, n = 12, Rank = "Species")

write.tidy(OTUTable, "Top12Species.txt")

#Plot
order <- OTUTable %>% filter(Species == "Bradyrhizobium sp. BTAi1") %>% arrange(RelativeAbundance)%>% extract2("Sample")
OTUTable <- OTUTable %>% mutate(Sample = factor(Sample, rev(order)))

Plot <- ggplot(OTUTable, aes(x = Sample, y = RelativeAbundance, fill = Species))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + facet_wrap(~Location,scales = "free_x")+ scale_y_continuous(expand = c(0,0))
Plot <- Plot + theme(axis.text.x = element_blank())
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + ylab(paste0("Relative Abundance (%)"))
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme(axis.title = element_text(size=16))
Plot <- Plot + theme(legend.title = element_text(size=12))
Plot <- Plot + theme(legend.text = element_text(size=10))
ggsave("top12species_by_sample_type.pdf",dpi=1086)
