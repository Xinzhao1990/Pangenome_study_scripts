library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)
library(reshape2)
source("AICc.PERMANOVA2.r")

tax <- read_tsv("kraken_species_abundance.tidy.txt")
tax <- as.data.frame(tax)
tax_cast <- dcast(tax, Sample ~ Species, value.var = "count", fill = 0)
colnames(tax_cast) <- c("Sample", paste0("Species",1:11784))

write_tsv(tax_cast, "kraken_species_abundance_rename.cast.txt")

#replace " " with "_" before loading the table
tax_cast <- read.table("kraken_species_abundance_pcoa_rename.cast.txt",header=T,row.names=1)
#distance <- distance(tax_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(tax_cast, method="bray")
distance <- vegdist(tax_cast, method="jaccard")

UniFracMatrix <- as.matrix(distance)
UniFrac <- melt(UniFracMatrix, value.name = "Distance")
names(UniFrac)[1:2] <- c("Sample1", "Sample2")

meta <- read_tsv("pan_metadata.txt") %>% select(Sample,Type)
names(meta)[1] <- c("Sample")

UniFrac <- merge(UniFrac, meta, by.x = "Sample1", by.y = "Sample", all.x = TRUE)
UniFrac <- merge(UniFrac, meta, by.x = "Sample2", by.y = "Sample", all.x = TRUE)
names(UniFrac)[4:5] <- c("Type1", "Type2")

UniFrac <- UniFrac[which( ! UniFrac$Sample1 == UniFrac$Sample2), ]

UniFrac$Comparison <- paste0(UniFrac$Type1, " vs. ", UniFrac$Type2)
UniFrac2 <- UniFrac %>% filter(Comparison %in% c("Skin vs. Skin", "Skin vs. Indoor surface", "Skin vs. Subway air", "Skin vs. Subway surface", "Skin vs. Urban public surface", "Skin vs. Pier surface"))
UniFrac2 <- UniFrac %>% filter(Comparison %in% c("Skin vs. Indoor surface", "Skin vs. Subway air", "Skin vs. Subway surface", "Skin vs. Urban public surface", "Skin vs. Pier surface"))

#UniFrac <- read_tsv("inter_type_bray_curtis_dissimilarity.txt")

library(viridis)
#Plot <- ggplot(UniFrac, aes(x = Comparison, y = Distance, fill = Comparison))
Plot <- ggplot(UniFrac2, aes(x = reorder(Comparison, Distance, median), y = Distance, fill = Comparison))
#Plot <- Plot + geom_boxplot()
Plot <- Plot + geom_violin(width=0.7) + geom_boxplot(width=0.1, color="#999999", alpha=0.2) + scale_fill_viridis(discrete = TRUE)
Plot <- Plot + theme_bw()
Plot <- Plot + scale_fill_brewer(palette = "Accent")
#Plot <- Plot + facet_wrap(~alpha,scales="free")
Plot <- Plot + theme(legend.position = "none")
Plot <- Plot + theme(axis.title = element_text(size=14, face="bold"))
Plot <- Plot + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
#Plot <- Plot + xlab(paste0("Pairwise comparison")) + ylab(paste0("Bray-Curtis dissimilarity"))
Plot <- Plot + xlab(paste0("Pairwise comparison")) + ylab(paste0("Jaccard distance"))

test <- wilcox.test(UniFrac$Distance~UniFrac$Comparison, data=UniFrac)
Wilcoxon rank sum test with continuity correction

data:  UniFrac$Distance by UniFrac$Comparison
W = 8707650, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0

test$p.value
[1] 0

library(pgirmess)
kruskalmc(UniFrac$Distance~UniFrac$Comparison, data=UniFrac)

#adonis <- adonis2(UniFracMatrix~Surface_material + Sample_type + Group + Location + Sample_type:Location,data=Meta)

PCoA <- cmdscale(UniFracMatrix, k = 2, eig = TRUE)
DF <- data.frame(Sample = row.names(PCoA$points), PCoA1 = PCoA$points[,1], PCoA2 = PCoA$points[,2], row.names = NULL)
Eigenvalues <- eigenvals(PCoA)
VarianceExplained <- Eigenvalues /sum(Eigenvalues)
VarianceExplained1 <- 100 * signif(VarianceExplained[1], 2)
VarianceExplained2 <- 100 * signif(VarianceExplained[2], 2)
PCoA <- merge(DF, Meta, by = "Sample", all.x = TRUE)

#Merge Dataframe with metadata
AllSamples <-data.frame(Sample = row.names(UniFracMatrix))
AllSamples <- merge(AllSamples, Meta, by = "Sample", all.x = TRUE)
#Check that UniFrac matrix rows match samples table (for ANOSIM)
sum(row.names(UniFracMatrix)==AllSamples$Sample) == length(AllSamples$Sample)

#Perform ANOSIM
BrayCurtis <- as.dist(distance)
ANOSIM <- anosim(BrayCurtis, grouping = AllSamples$Location)
#See ANOSIM options
ls(ANOSIM)

#Load ANOSIM statistics and signif
ANOSIM$statistic
ANOSIM$signif

#Perform Adonis(permanova)
Adonis <- adonis2(BrayCurtis~Surface_material + Sample_type + Group + Location + Sampling_date + Air_temperature + Air_humidity+ Location:Sample_type+Location:Surface_material+Location:Air_temperature+Location:Air_humidity+Location:Sampling_date+Location:Group+ Sample_type:Surface_material+Sample_type:Air_temperature+Sample_type:Air_humidity+Sample_type:Sampling_date+Sample_type:Group+ Surface_material:Air_temperature+Surface_material:Air_humidity+Surface_material:Sampling_date+Surface_material:Group+ Air_temperature:Air_humidity+Air_temperature:Sampling_date+Air_temperature:Group+ Air_humidity:Sampling_date+Air_humidity:Group+ Sampling_date:Group, data=AllSamples)


#best model for nested permanova
Adonis <- adonis2(BrayCurtis~Surface_material + Sample_type + Group + Location + Sample_type:Location,data=AllSamples)


AICc.PERMANOVA2(Adonis)$AIC
[1] -484.7012

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = BrayCurtis ~ Surface_material + Sample_type + Group + Location + Sample_type:Location, data = AllSamples)
Df SumOfSqs      R2        F Pr(>F)
Surface_material       1    6.139 0.17157 117.3879  0.001 ***
Sample_type            2    2.250 0.06287  21.5079  0.001 ***
Group                  1    2.906 0.08123  55.5742  0.001 ***
Location               7    4.697 0.13127  12.8302  0.001 ***
Sample_type:Location  24   12.519 0.34990   9.9747  0.001 ***
Residual             139    7.269 0.20316
Total                174   35.780 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


dbrda <- dbrda(BrayCurtis ~ Sample_type+Location+Surface_material, data = AllSamples)
anova(dbrda,by="margin")


Adonis <- adonis(BrayCurtis~Location, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Location, data = AllSamples, permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Location    8     7.715 0.96444  5.7046 0.21564  0.001 ***
Residuals 166    28.064 0.16906         0.78436
Total     174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Sample_type, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Sample_type, data = AllSamples,      permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sample_type   3     8.388 2.79615  17.456 0.23445  0.001 ***
Residuals   171    27.391 0.16018         0.76555
Total       174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Surface_material, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Surface_material, data = AllSamples,      permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Surface_material   1     6.139  6.1389   35.83 0.17157  0.001 ***
Residuals        173    29.641  0.1713         0.82843
Total            174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Group, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Group, data = AllSamples, permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Group       1     2.892  2.8918  15.211 0.08082  0.001 ***
Residuals 173    32.888  0.1901         0.91918
Total     174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Specify X and Y axis
xlab <- paste0("PCoA1 (Variance Explained: 31%)")
ylab <- paste0("PCoA2 (Variance Explained: 12%)")

#Draw plot by sample types
library(ggthemes)

Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, colour=Group))
#Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Group))
Plot <- Plot + geom_point(size=3)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Set1")
Plot <- Plot + theme(axis.title=element_text(size=18, face="bold"))
Plot <- Plot + theme(axis.text=element_text(size=14, face="bold"))
ggsave("bray_curits_pcoa_by_north_south.pdf")
ggsave("bray_curits_pcoa_by_north_south_no_ellipse.pdf")


Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, colour=Location, shape=Sample_type))
#Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Sample_type))
Plot <- Plot + geom_point(size=2) + aes(fill = Location)
Plot <- Plot + facet_wrap(~Location)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
#Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Set1")
Plot <- Plot + theme(axis.title=element_text(size=12, face="bold"))
Plot <- Plot + theme(axis.text=element_text(size=10))
ggsave("bray_curits_pcoa_by_location.pdf")
ggsave("bray_curits_pcoa_by_location_facet.pdf")


Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, colour=Sample_type))
Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Sample_type))
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Set3")
Plot <- Plot + theme(axis.title=element_text(size=18, face="bold"))
Plot <- Plot + theme(axis.text=element_text(size=14, face="bold"))
ggsave("bray_curits_pcoa_by_sample_type.pdf")

Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, color=Surface_material, shape=Sample_type))
#Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Surface_material))
Plot <- Plot + geom_point(size=3)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title=element_text(size=18, face="bold"))
Plot <- Plot + theme(axis.text=element_text(size=14, face="bold"))
ggsave("bray_curits_pcoa_by_material_type.pdf")
