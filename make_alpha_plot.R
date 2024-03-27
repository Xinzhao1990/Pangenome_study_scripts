library(readr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

#make the boxplot for richness
table <- read_tsv("pan_alpha_diversity_rarefied.txt")
table <- table %>%
gather(alpha, value, 2:(ncol(table))) %>%
filter(value > 0)

meta <- read_tsv("pan_metadata.txt")
table <- left_join(table,meta) %>% filter(!alpha == "chao1_richness")
table$Type <- factor(table$Type, levels = c("Skin","Indoor surface","Subway air", "Subway surface",  "Urban public surface","Pier surface"))
#Plot <- ggplot(table, aes(x = reorder(Sample, value, median), y = richness, fill = Type))
Plot <- ggplot(table, aes(Type,value, fill = Type))
Plot <- Plot + geom_boxplot()
Plot <- Plot + theme_bw()
Plot <- Plot + facet_wrap(~alpha,scales="free")
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme(axis.title = element_text(size=14))
#Plot <- Plot + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
Plot <- Plot + theme(axis.text.x = element_blank())
Plot <- Plot + theme(axis.ticks.x=element_blank())
Plot <- Plot + ylab(paste0("Value")) + xlab(paste0("Habitat"))
ggsave("pan_alpha_diversity_rarefied_by_sample_type.pdf", dpi=1086)


table <- read_tsv("pan_samples_kraken_subsamp_alpha_diversity.txt")
table$Type <- factor(table$Type, levels = c("Skin","Indoor surface","Subway air", "Subway surface", "Pier surface", "Urban public surface"))
Plot <- ggplot(table, aes(x = Type, y = value, fill = alpha))
Plot <- Plot + geom_boxplot()
#Plot <- Plot + geom_violin(width=0.7) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + scale_fill_viridis(discrete = TRUE)
Plot <- Plot + theme_bw()
#Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + facet_wrap(~alpha,scales="free")
Plot <- Plot + theme(axis.title = element_text(size=14, face="bold"))
Plot <- Plot + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
Plot <- Plot + ylab(paste0("Richness")) + xlab(paste0("Sample Type"))
ggsave("pier_samples_kraken_subsamp_richness_boxplot_sample_type_location_violin.pdf", dpi=1086)

table <- read_tsv("pier_samples_kraken_subsamp_richness.txt")
Plot <- ggplot(table, aes(x = Surface_material, y = richness, fill = Surface_material))
Plot <- Plot + geom_boxplot()
Plot <- Plot + theme_bw()
Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title = element_text(size=14, face="bold"))
#Plot <- Plot + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
Plot <- Plot + ylab(paste0("Richness")) + xlab(paste0("Surface Material"))
ggsave("pier_samples_kraken_subsamp_richness_boxplot_material.pdf", dpi=1086)

table <- read_tsv("pier_samples_kraken_subsamp_richness.txt")
Plot <- ggplot(table, aes(x = Group, y = richness, fill = Group))
Plot <- Plot + geom_boxplot()
Plot <- Plot + theme_bw()
Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title = element_text(size=14, face="bold"))
#Plot <- Plot + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
Plot <- Plot + ylab(paste0("Richness")) + xlab(paste0("Group"))
ggsave("pier_samples_kraken_subsamp_richness_boxplot_north_south.pdf", dpi=1086)


library(stats)
#Statistics
stats <- aov(richness~Location, data=table)
summary(stats)
Df  Sum Sq Mean Sq F value   Pr(>F)
Location      8 3285174  410647   8.951 6.06e-10 ***
Residuals   146 6698123   45878
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

stats <- aov(richness~Sample_type, data=table)
Df  Sum Sq Mean Sq F value   Pr(>F)
Sample_type   3 1759354  586451   10.77 1.87e-06 ***
Residuals   151 8223943   54463
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

stats <- aov(richness~Surface_material, data=table)
Df  Sum Sq Mean Sq F value   Pr(>F)
Surface_material   1 1155900 1155900   20.04 1.48e-05 ***
Residuals        153 8827397   57695
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

stats <- aov(richness~Setting, data=table)
> summary(stats)
Df  Sum Sq Mean Sq F value Pr(>F)
Setting       1  133891  133891    2.08  0.151
Residuals   153 9849406   64375





