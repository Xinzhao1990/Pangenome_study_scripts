library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)
library(reshape2)
source("AICc.PERMANOVA2.r")

tax <- read_tsv("kraken_report_all_species_w_unclassified.tidy.txt")
tax <- tax %>% filter(!Species == "unclassified")
tax <- as.data.frame(tax)

tax_cast <- dcast(tax, Sample ~ Species, value.var = "Count", fill = 0)
colnames(tax_cast) <- c("Sample", paste0("Species",1:8250))

write_tsv(tax_cast, "kraken_species_abundance_pcoa_rename.cast.txt")

#replace " " with "_" before loading the table
tax_cast <- read.table("kraken_species_abundance_pcoa_rename.cast.txt",header=T,row.names=1)
#distance <- distance(tax_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(tax_cast, method="bray")
distance <- vegdist(tax_cast, method="jaccard")

UniFracMatrix <- as.matrix(distance)
write.table(UniFracMatrix,"kraken_species_bray_cutis_matrix.txt", sep="\t", col.names=TRUE, row.names=TRUE)

Meta <- read_tsv("pan_metadata.txt")

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
Adonis <- adonis2(BrayCurtis~Type,data=AllSamples)
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = BrayCurtis ~ Type, data = AllSamples)
Df SumOfSqs      R2      F Pr(>F)    
Type       5    37.52 0.21247 39.498  0.001 ***
  Residual 732   139.07 0.78753                  
Total    737   176.59 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(pairwiseAdonis)
# Pairwise PERMANOVA
pairwise_results <- pairwise.adonis(BrayCurtis, AllSamples$Type)
print(pairwise_results)


#best model for nested permanova
Adonis <- adonis2(BrayCurtis~Surface_material + Sample_type + Group + Location + Sample_type:Location,data=AllSamples)


AICc.PERMANOVA2(Adonis)$AIC
[1] -484.7012

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

pairwise.perm.manova(BrayCurtis~AllSamples$Type)

# calculate multivariate dispersions
mod <- betadisper(BrayCurtis, AllSamples$Type)
mod

## Perform test
anova(mod)
Analysis of Variance Table

Response: Distances
Df Sum Sq Mean Sq F value    Pr(>F)    
Groups      5 4.5995 0.91989  89.514 < 2.2e-16 ***
  Residuals 732 7.5224 0.01028                      
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pmod <- permutest(mod, permutations = 99, pairwise = TRUE)

# perform a permutation-based test.
adonis2(dist(mod$distances)~AllSamples$Type)
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist(mod$distances) ~ AllSamples$Type)
Df SumOfSqs      R2      F Pr(>F)    
AllSamples$Type   5   4.5995 0.37944 89.514  0.001 ***
  Residual        732   7.5224 0.62056                  
Total           737  12.1219 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 99

Response: Distances
Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
Groups      5 4.5995 0.91989 89.514     99   0.01 **
  Residuals 732 7.5224 0.01028                        
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Pairwise comparisons:
  (Observed p-value below diagonal, permuted p-value above diagonal)
Indoor surface        Skin Subway air Subway surface Urban public surface
Indoor surface                        3.0000e-02 1.5000e-01 1.0000e-02     1.0000e-02                 0.06
Pier surface             3.1405e-02              6.9000e-01 1.0000e-02     1.0000e-02                 0.01
Skin                     1.1714e-01   6.8329e-01            1.0000e-02     1.0000e-02                 0.01
Subway air               1.7325e-40   9.1188e-50 1.5918e-40                1.0000e-02                 0.01
Subway surface           2.2194e-13   5.3836e-19 4.6924e-15 3.2844e-08                                0.01
Urban public surface     3.5625e-02   5.8414e-05 9.7750e-04 1.0996e-27     3.1874e-07    

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
#bray-curtis
xlab <- paste0("PCoA1 (Variance explained: 20%)")
ylab <- paste0("PCoA2 (Variance explained: 13%)")
#jaccard
xlab <- paste0("PCoA1 (Variance explained: 13%)")
ylab <- paste0("PCoA2 (Variance explained: 8.1%)")

#Draw plot by sample types
library(ggthemes)

PCoA$Type <- factor(PCoA$Type, levels = c("Skin","Indoor surface","Subway air","Subway surface", "Urban public surface","Pier surface"))
Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, color=Type,shape=Type))
Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Type))
#Plot <- Plot + geom_point(size=2)
#Plot <- Plot + facet_wrap(~Sample_type)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Paired")
Plot <- Plot + theme(axis.title=element_text(size=12))
Plot <- Plot + theme(axis.text=element_text(size=10))
ggsave("bray_curtis_pcoa_by_sample_type_ellipse.pdf")
ggsave("jaccard_pcoa_by_sample_type_ellipse.pdf")

Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, color=Category))
Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Category))
Plot <- Plot + geom_point(size=2)
#Plot <- Plot + facet_wrap(~Sample_type)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Set2")
Plot <- Plot + theme(axis.title=element_text(size=12))
Plot <- Plot + theme(axis.text=element_text(size=10))
ggsave("bray_curits_pcoa_by_project_ellipse.pdf")
ggsave("jaccard_pcoa_by_project_ellipse.pdf")


Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, color=Surface_material,shape=Surface_material))
Plot <- Plot + geom_point(size=2) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Surface_material))
Plot <- Plot + geom_point(size=2)
#Plot <- Plot + facet_wrap(~Sample_type)
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_color_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title=element_text(size=12))
Plot <- Plot + theme(axis.text=element_text(size=10))
ggsave("bray_curits_pcoa_by_material_ellipse.pdf")
ggsave("jaccard_pcoa_by_material_ellipse.pdf")

Plot2 <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2, color=Surface_material))
Plot2 <- Plot2 + geom_point(size=2,aes(shape=Sample_type)) + stat_ellipse(geom = "polygon",alpha = 0.05, aes(fill = Surface_material))
#Plot <- Plot + geom_point(size=3)
Plot2 <- Plot2 + xlab(xlab)
Plot2 <- Plot2 + ylab(ylab)
Plot2 <- Plot2 + theme_bw()
Plot2 <- Plot2 + scale_color_brewer(palette = "Dark2")
Plot2 <- Plot2 + theme(axis.title=element_text(size=12))
Plot2 <- Plot2 + theme(axis.text=element_text(size=12))
Plot2 <- Plot2 + ggtitle("Unrarefied dataset")+theme(plot.title = element_text(hjust = 0.5))
ggsave("bray_curits_pcoa_by_material_type_ellipse.pdf")

library(ggpubr)
library(gridExtra)
Figure <- grid.arrange(Plot2, Plot, ncol=2)
ggsave("pier_samples_kraken_subsamp_shannon_boxplot_combined_by_sample_type_location.pdf",Figure)
