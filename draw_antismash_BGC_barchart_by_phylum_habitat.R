library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(wilkoxmisc)
library(reshape2)
library(RColorBrewer)

phylum <- read_tsv("BGC_percentage_by_phylum_gtdbtk.txt") %>%
select(BGC,percentage,Phylum_count,Category)
names(phylum)[3] <- c("Group")
location <- read_tsv("BGC_percentage_by_location_gtdbtk.txt")%>%
select(BGC,percentage,Location_count,Category)
names(location)[3] <- c("Group")
merge <- rbind(phylum,location)

type <- read_tsv("BGC_percentage_by_type.txt")%>%
select(BGC,percentage,Type_count,Category)
names(type)[3] <- c("Group")
merge <- rbind(merge,type)

merge2 <- read_tsv("BGC_percentage_by_phylum_habitat_combined.txt")

colorCount = length(unique(merge2$BGC))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

cbPalette2 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c', '#fb9a99','#e31a1c', '#fdbf6f','#ff7f00', '#cab2d6', '#6a3d9a','#ffff99', '#b15928','#d9d9d9')


merge2$Category <- factor(merge2$Category,levels=c("Phylum","Type"))
Plot <- ggplot(merge2, aes(y = Group, x = percentage, fill=BGC))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + theme_classic()
Plot <- Plot + facet_wrap(Category~.,scales = "free_y")
Plot <- Plot + scale_fill_manual(values= getPalette(colorCount))
#Plot <- Plot + scale_fill_manual(values = cbPalette2)
#Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + theme(legend.position = "bottom")
#Plot <- Plot + theme(axis.text.y = element_text(face="italic"))
Plot <- Plot + xlab(paste0("Percentage (%)"))


Plot <- ggplot(merge, aes(y = reorder(BGC, number), x = number, fill=BGC))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + scale_fill_manual(values= getPalette(colorCount))
#Plot <- Plot + coord_flip()
Plot <- Plot + facet_wrap(~Category)
Plot <- Plot + theme_classic()
#Plot <- Plot + theme(axis.text.x = element_text(angle=90,hjust=1,size=6))
#Plot <- Plot + theme(axis.text.y = element_text(angle=90,hjust=1,size=6))
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + xlab(paste0("BGC types")) + ylab(paste0("Count"))
Plot <- Plot + theme(legend.position = "none")

ggsave("BGC_count_by_type_barchart_gtdbtk_new.pdf")
#make the barchart by phylum (average BGC count per phylum)
table <- read_tsv("antismash_BGC_count_per_bin.txt")
bgc <- read_tsv("BGC_regroup.txt")
merge <- left_join(table,bgc) %>%
select(BGC,Bin,number) %>%
group_by(Bin,BGC) %>%
mutate(number=sum(number)) %>%
ungroup() %>%
unique()

genome <- read_tsv("gtdbtk.bac120.summary.txt")
names(genome)[1] <- c("Bin")
genome$Bin <- gsub("bin.","bin", genome$Bin)

merge2 <- left_join(merge,genome)%>%
select(BGC,number,Phylum) %>%
group_by(BGC,Phylum) %>%
mutate(number=sum(number)) %>%
ungroup() %>%
unique() %>%
group_by(Phylum) %>%
mutate(sum=sum(number))%>%
ungroup() %>%
mutate(percentage=number*100/sum)

merge2 <- merge2 %>%
mutate(Phylum = gsub("p__", "", merge2$Phylum))

write_tsv(merge2,"BGC_percentage_by_phylum_gtdbtk.txt")

merge2 <- read_tsv("BGC_percentage_by_phylum_gtdbtk.txt")

colorCount = length(unique(merge2$BGC))
getPalette = colorRampPalette(brewer.pal(12, "Set1"))

Plot <- ggplot(merge2, aes(y = Phylum_count, x = percentage, fill=BGC))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + theme_classic()
Plot <- Plot + scale_fill_manual(values= getPalette(colorCount))
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + xlab(paste0("Percentage (%)")) + ylab(paste0("Phylum"))

#make the barchart by location (average BGC count per location)
table <- read_tsv("antismash_BGC_count_per_bin.txt")
bgc <- read_tsv("BGC_regroup.txt")
table <- left_join(table,bgc) %>%
select(BGC,Bin,number) %>%
group_by(Bin,BGC) %>%
mutate(number=sum(number)) %>%
ungroup() %>%
unique()

table$Bin <- gsub("_"," ",table$Bin)
table <- table %>%
    mutate(Sample=word(table$Bin,1))
meta <- read_tsv("pier_metadata_complete.txt") %>%
    select(Sample,Location, Sample_type,Surface_material)

merge <- left_join(table,meta) %>%
select(BGC,number,Location) %>%
group_by(BGC,Location) %>%
mutate(number=sum(number)) %>%
ungroup() %>%
unique() %>%
group_by(Location) %>%
mutate(sum=sum(number))%>%
ungroup() %>%
mutate(percentage=number*100/sum)

merge2 <- read_tsv("BGC_percentage_by_location_gtdbtk.txt")


colorCount = length(unique(merge$BGC))
getPalette = colorRampPalette(brewer.pal(12, "Set1"))

Plot <- ggplot(merge2, aes(y = Location_count, x = percentage, fill=BGC))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + scale_fill_manual(values= getPalette(colorCount))
#Plot <- Plot + coord_flip()
#Plot <- Plot + theme_minimal()
#Plot <- Plot + theme(axis.text.x = element_text(angle=90,hjust=1,size=6))
#Plot <- Plot + theme(axis.text.y = element_text(angle=90,hjust=1,size=6))
Plot <- Plot + theme(axis.ticks.x = element_blank())
Plot <- Plot + xlab(paste0("Percentage (%)")) + ylab(paste0("Location"))





