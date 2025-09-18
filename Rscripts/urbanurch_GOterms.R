#TopGO analysis with ubanurchin outliers
#Madison Armstrong
#8/13/25

setwd("~Desktop/urbanurchins")

#TopGO script
#Madison Armstrong 
#Last Modified 7/7/2025

#Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(fuzzyjoin)

#read in files####
#https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000002235.5/ base datasheet downloaded from here
ncbi_all<-read_tsv("genomicdata/ncbi_dataset.tsv")

ncbi_all<-ncbi_all %>% 
  select(c("Accession","Begin","End", "Name", "Gene ID", "Gene Type","Transcripts accession", "Transcript name")) #organize

ncbi_all<-ncbi_all %>% 
  separate(col=Accession, into= c("Chunk", "Version"), sep="\\.")

#shared SNPs--pulled from outlier script
overlapping_snps<-read.csv("genomicdata/overlapping_snps_all.csv")

#separate accession and location for later merging with larger ncbi dataset
overlapping_snps<-overlapping_snps %>% 
  separate(col=overlaps, into= c("Accession", "Location"), sep=".1_")


###outliers for each regions####
sigVic<-read.csv("genomicdata/Vic_outliers.csv") #17 outliers
sigVic<-sigVic %>% 
  separate(col=LocusName, into= c("Accession", "Location"), sep=".1_")

sigLA<-read.csv("genomicdata/LA_outliers.csv") #2271 outliers
sigLA<-sigLA %>% 
  separate(col=LocusName, into= c("Accession", "Location"), sep=".1_")

sigSD<-read.csv("genomicdata/SD_outliers.csv") #749 outliers
sigSD<-sigSD %>% 
  separate(col=LocusName, into= c("Accession", "Location"), sep=".1_")

###first identify GO terms for shared overlapping SNPs, but need to fuzzy join since I am looking for a range
merged_overlap <- overlapping_snps %>%
  inner_join(ncbi_all, by = c("Accession" = "Chunk")) %>%
  filter(Location >= Begin & Location <= End)

#keep only distinct rows
merged_overlap <- merged_overlap %>%
  distinct(Location, .keep_all = TRUE)
write.csv(merged_overlap, "genomicdata/overlap_GOTerms.csv")

###for each region!####
#Vic
merged_Vic <- sigVic %>%
  inner_join(ncbi_all, by = c("Accession" = "Chunk")) %>%
  filter(Location >= Begin & Location <= End)
#keep only distinct rows
merged_Vic <- merged_Vic %>%
  distinct(Location, .keep_all = TRUE) #16 GO Terms
write.csv(merged_Vic, "genomicdata/Vic_GOTerms.csv")

#LA
merged_LA <- sigLA %>%
  inner_join(ncbi_all, by = c("Accession" = "Chunk")) %>%
  filter(Location >= Begin & Location <= End)
#keep only distinct rows
merged_LA <- merged_LA %>%
  distinct(Location, .keep_all = TRUE) #2034 GO Terms
write.csv(merged_LA, "genomicdata/LA_GOTerms.csv")

#SD
merged_SD <- sigSD %>%
  inner_join(ncbi_all, by = c("Accession" = "Chunk")) %>%
  filter(Location >= Begin & Location <= End)
#keep only distinct rows
merged_SD <- merged_SD %>%
  distinct(Location, .keep_all = TRUE) #670 GO Terms
write.csv(merged_SD, "genomicdata/SD_GOTerms.csv")

