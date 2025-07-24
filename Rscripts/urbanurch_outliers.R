#outlier analyses
#Madison Armstrong
#2/14/25

setwd("~Desktop/urbanurchins")
# #install via 
#if (!("devtools" %in% installed.packages())){install.packages("devtools")}
# 
# ### Install required package qvalue if not installed
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
# 
#BiocManager::install("pcadapt")
#install.packages("pcadapt")
# 
# ### Install OutFLANK if not installed
#devtools::install_github("whitlock/OutFLANK")
#if (!("vcfR" %in% installed.packages())){install.packages("vcfR")}

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "white", fill = NA, size = 1)) 
}

library(devtools)
library(tidyverse)
library(OutFLANK)
library(vcfR)
library(pcadapt)
library(dplyr)
library(patchwork)
library(BEDMatrix)
library(bigsnpr)
library(data.table)
library(eulerr)

#remotes::install_github("kaustubhad/fastman", ref = "main", force = TRUE) #had to download this way
library(fastman)

####Prepping bed, fam and bim files####
#clean and convert vcf to bed file in terminal
# Path to your PLINK files (without extension)
#snp_readBed("genomicdata/7.urbanurchfull.new.bed")
# Load as bigSNP object
obj <- snp_attach("genomicdata/7.urbanurchfull.new.rds")
G <- obj$genotypes      # Genotype matrix
CHR <- obj$map$chromosome  # Chromosome info
SAMPLES <- obj$fam      # Sample info

# Convert bigSNP object to a regular matrix
G_matrix <- G[]

# Merge data
merged_data <- cbind(SAMPLES, G_matrix)
colnames(merged_data)[-(1:ncol(SAMPLES))] <- CHR  # Assign chromosome labels to SNPs
head(merged_data[,1:10])
#drop unnecessary columns
merged_dataset <- subset(merged_data, select = -c(3:6) )
head(merged_dataset[,1:10]) #now we can work with this!


####PCADAPT####
##pcadapt from a bed file
genotypes <- read.pcadapt("genomicdata/7.urbanurchfull.new.bed", type= "bed")  
pcadapt <- pcadapt(input = genotypes, K =6) #since later pc axes are interesting!!
pcadapt.comp <- pcadapt(input = genotypes, K =6, method="componentwise")

##Variance explained per axis
EV <- pcadapt$singular.values^2
EV

#####visualization plots####
plot(pcadapt, option = "screeplot") 
#manhattan plot
plot(pcadapt, option="manhattan")
#qqplot
plot(pcadapt, option="qqplot")
#histogram
hist(pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#dist
plot(pcadapt, option = "stat.distribution")

#####Read in metadata#######

#read in .bim file
bim_data <- read.table("genomicdata/7.urbanurchfull.new.bim", header = FALSE,sep = "\t")[,c(1,4)]
names(bim_data) <- c("Chr","Pos")
bim_data$snp_names <- paste(bim_data$Chr,"_",bim_data$Pos,sep="")
bim_data$Chr <- as.character(bim_data$Chr)

inds <- read.table("genomicdata/7.urbanurchfull.new.fixed.fam") #fixed file in text editor so samples don't have a space
meta <- read.csv("metadata/urbanurchins_metadata_thinned.csv")
meta$Sample_ID <- gsub(" ", "", meta$Sample_ID) # some samples have a space, making them not match later
meta.order <- meta[match(inds$V1,meta$Sample_ID),] #183 samples now woo


##interaction of region and urban/nonurban boxplots
PC3 <- ggplot(data=meta.order,aes(x=meta.order$region,y=pcadapt$scores[,3],fill=meta.order$Dev))  + xlab("Region")+ ylab("Pcadapt Scores (PC3)")+ geom_boxplot()+theme_box()+ scale_fill_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))
PC5 <- ggplot(data=meta.order,aes(x=meta.order$region,y=pcadapt$scores[,5],fill=meta.order$Dev)) + xlab("Region")+ ylab("Pcadapt Scores (PC5)") + geom_boxplot()+theme_box()+ scale_fill_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))
PC6 <- ggplot(data=meta.order,aes(x=meta.order$region,y=pcadapt$scores[,6],fill=meta.order$Dev)) + xlab("Region")+ ylab("Pcadapt Scores (PC6)") + geom_boxplot() +theme_box() + scale_fill_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))
pcadapt_int<-
  PC3+PC5+PC6 +plot_layout(guides="collect",nrow=1)

ggsave("PCadapt_interaction.png", pcadapt_int, width=40, height=15, units = "cm")

##Make full pcadapt frame
pcframe <- data.frame(bim_data,PC.full=p.adjust(pcadapt$pvalues,method="fdr"),
                      PC1=p.adjust(pcadapt.comp$pvalues[,1],method="fdr"),
                      PC2=p.adjust(pcadapt.comp$pvalues[,2],method="fdr"),
                      PC3=p.adjust(pcadapt.comp$pvalues[,3],method="fdr"),
                      PC4=p.adjust(pcadapt.comp$pvalues[,4],method="fdr"),
                      PC5=p.adjust(pcadapt.comp$pvalues[,5],method="fdr"),
                      PC6=p.adjust(pcadapt.comp$pvalues[,6],method="fdr"))
manhattan_plot(pcadapt.comp,chr.info=as.numeric(bim_data$Chr),snp.info=bim_data$snp_names,K=3)
pcframe.filt <- pcframe %>% filter(!is.na(PC.full))
dim(pcframe.filt) #V1 scaffold, V4 position

#trim short contigs
contig_counts <- table(pcframe.filt$Chr) #helpful visualization

min_snps <- 400  #all others in the 1000s

pcframe_bigcontig <- pcframe.filt %>%
  group_by(Chr) %>%
  filter(n() >= min_snps) %>%
  ungroup()
bigcontig_counts <- table(pcframe_bigcontig$) #helpful visualization
#KT<-write.csv(bigcontig_counts, 'bigcontig_counts.csv')


png("full.manhattan_plot.png", width = 1000, height = 400, res = 150)
fastman(pcframe_bigcontig,chr="V1",bp="V4", p="PC.full")
dev.off() #closes graphics device and writes file

png("manhattan_plotPC3.png", width = 1000, height = 600, res = 150)
fastman(pcframe_bigcontig,chr="V1",bp="V4",p="PC3")
dev.off() #closes graphics device and writes file

png("manhattan_plotPC5.png", width = 1000, height = 600, res = 150)
fastman(pcframe_bigcontig,chr="V1",bp="V4",p="PC5", suggestiveline=0.1)
dev.off() #closes graphics device and writes file

png("manhattan_plotPC6.png", width = 1000, height = 600, res = 150)
fastman(pcframe_bigcontig,chr="V1",bp="V4",p="PC6", suggestiveline=0.1)
dev.off() #closes graphics device and writes file


#stats for PC axes
pcdata<-cbind(meta.order, pcadapt$scores)
head(pcdata)
#PC1 not sig
#PC2
pc2<-lm(pcdata$'2'~pcdata$region+pcdata$Dev+pcdata$Depth..m.)
summary(pc2) #slight for region (Vic) and depth (0.02)
#PC3
pc3<-lm(pcdata$'3'~pcdata$region+pcdata$Dev+pcdata$Depth..m.)
summary(pc3) #sig by urban/nonurban!! (p=0.004) and depth (p=0.004)
#PC4
pc4<-lm(pcdata$'4'~pcdata$region+pcdata$Dev+pcdata$Depth..m.)
summary(pc4) #slight for SD & vic sig for depth (p<0.001) 
#PC5
pc5<-lm(pcdata$'5'~pcdata$region+pcdata$Dev+pcdata$Depth..m.)
summary(pc5) #sig by urban/nonurban (p<0.001), region for Victoria(makes sense!)
#PC6
pc6<-lm(pcdata$'6'~pcdata$region+pcdata$Dev+pcdata$Depth..m.)
summary(pc6) #sig by urban/nonurban, depth & region for SD

#####PC plots for different comparisons####
scores <- as.data.frame(pcadapt$scores)
scores$pop <- factor(meta.order$Dev)
colnames(scores)[1:6] <- paste0("PC", 1:6)


#dev + stat_ellipse(data=meta.order$Dev)
pl <- ggplot(scores, aes(x = PC1, y = PC2, color = pop, fill = pop)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.15, color = NA) + stat_ellipse(geom = "path", size = 1) +
  scale_color_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  scale_fill_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  theme_box()

p2 <- ggplot(scores, aes(x = PC3, y = PC4, color = pop, fill = pop)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.15, color = NA) + stat_ellipse(geom = "path", size = 1) +
  scale_color_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  scale_fill_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  theme_box()

p3 <- ggplot(scores, aes(x = PC5, y = PC6, color = pop, fill = pop)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.15, color = NA) + stat_ellipse(geom = "path", size = 1) +
  scale_color_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  scale_fill_manual(values = c("nonurban" = "#99d1ec", "urban" = "#665d4b")) +
  theme_box()

# Combine the plots in a 1x3 grid
pca_urb<-(pl+p2+p3 + plot_layout(ncol = 3, guides="collect"))
ggsave("pcadapt_urb.png", pca_urb, width=35, height=15, units = "cm")

#region
scores$region <- meta.order$region
p11<-ggplot(scores, aes(x = PC1, y = PC2, color = region, fill = region)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) +stat_ellipse(geom = "path", size = 1)+ 
  scale_color_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2")) +
  scale_fill_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2"))+theme_box()

p22<-ggplot(scores, aes(x = PC3, y = PC4, color = region, fill = region)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) + stat_ellipse(geom = "path", size = 1) +
  scale_color_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2"))+
  scale_fill_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2")) +theme_box()
p33<-ggplot(scores, aes(x = PC5, y = PC6, color = region, fill = region)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) + stat_ellipse(geom = "path", size = 1)+ 
  scale_color_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2")) +
  scale_fill_manual(values = c("Vic"="orange" , "LA" = "purple2",  "SD" = "orchid2")) +theme_box()
# Combine the plots in a 1x3 grid
pca_region<-(p11/p22/p33)
ggsave("pcadapt_region.png", pca_region, width=30, height=20, units = "cm")


#tidal
scores$tidal <- meta.order$Depth..m.
a<-ggplot(scores, aes(x = PC1, y = PC2, color = tidal, fill = tidal)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) +
  stat_ellipse(geom = "path", size = 1)+
  scale_color_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) + 
    scale_fill_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) +theme_box()

b<-ggplot(scores, aes(x = PC3, y = PC4, color = tidal, fill = tidal)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) +
  stat_ellipse(geom = "path", size = 1)+
  scale_color_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) + scale_fill_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) +theme_box()

c<-ggplot(scores, aes(x = PC5, y = PC6, color = tidal, fill = tidal)) +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, color = NA) +
  stat_ellipse(geom = "path", size = 1) +
  scale_color_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) +scale_fill_manual(values=c("intertidal"="tan", "subtidal"="darkblue")) +theme_box()
# Combine the plots in a 1x3 grid
pca_tidal<-(a/b/c)
ggsave("pcadapt_tidal.png", pca_tidal, width=30, height=20, units = "cm")


###code doesn't rerun from here #####
#####identifying outliers & making csv 
#outlier line, plot manhattan plots-pvalues
alpha <- 0.1  # Desired FDR level
adjusted_pvalues <- p.adjust(pcadapt_snp_info$P_value, method = "fdr")  # Adjust p-values
#identify outliers
pcadapt_outliers <- which(adjusted_pvalues < alpha) #still in order numerically
length(pcadapt_outliers) #191!!!

pcadapt_outliers_info <- data.frame(
  LocusName = pcadapt_snp_info$SNP_Name[pcadapt_outliers],
  P_value = pcadapt_snp_info$P_value[pcadapt_outliers],
  Adjusted_P_value = adjusted_pvalues[pcadapt_outliers]
)
write.csv(pcadapt_outliers_info, "pcadapt_outliers_info.csv")
pcadapt_names<-pcadapt_outliers_info$LocusName
write.csv(pcadapt_names, "pcadapt_names.csv")


#for later manhattan plot
threshold <- -log10(max(x$pvalues[adjusted_pvalues < alpha], na.rm = TRUE))  # Get -log10 threshold

#####Association of PCs with outliers:
snp_pc <- get.pc(x, pcadapt_outliers)
# Identify which outliers are most associated with PC2 or PC3
pc2_outliers <- pcadapt_outliers[snp_pc$PC == 2] #22 SNPs
pc3_outliers <- pcadapt_outliers[snp_pc$PC == 3] #15 SNPs
pc4_outliers <- pcadapt_outliers[snp_pc$PC == 4] #28 SNPs
pc5_outliers <- pcadapt_outliers[snp_pc$PC == 5] #10 SNPs
pc6_outliers <- pcadapt_outliers[snp_pc$PC == 6] #20 SNPs
pc7_outliers <- pcadapt_outliers[snp_pc$PC == 7] #21 SNPs
# Create a color vector based on whether the SNP is in combined_outliers
colors <- rep("grey", length(x$pvalues))  # Initialize all colors to grey
if (length(pcadapt_outliers) > 0) {
  colors[pc2_outliers] <- "pink" 
  colors[pc3_outliers] <-"orange" 
  colors[pc4_outliers] <-"darkred" 
  colors[pc5_outliers] <-"purple" 
  colors[pc6_outliers] <-"blue" 
  colors[pc7_outliers] <-"darkgreen" 
}
#####manhattan plot#
# Plot all SNPs but PC2 and PC3 ones in grey!
plot(-log10(x$pvalues), pch = 19, cex = 0.3, col = colors, 
     ylab = "-log10 P-value", xlab = "SNPs")
# add line
abline(h = threshold, col = "red", lty = 2)
legend("bottomright", legend=c("PC2", "PC3", "PC4", "PC5", "PC6", "PC7"), col=c("pink", "orange", "darkred", "purple", "blue","darkgreen"), pch=19,horiz=TRUE)




###OUTFLANK#####
# read in bedfile
#read in with bed2matrix in pcadapt converts bed file to matrix
path_to_file <- "genomicdata/7.urbanurchfull.new.bed"
genos<-bed2matrix(path_to_file) #.fam file must exist for this step to run
#take a look at it
tail(genos[,1:10]) #some NAs!
#genotypes need to be 0, 1 or 2 so all okay here!
genos[is.na(genos)]<- 9 #change NA's to 9 for outflank
dim(genos) # 183 219773 looks good!

#read in metadata
sites.thinned<-read.csv('metadata/urbanurchins_metadata_thinned.csv', header=TRUE, sep=',')
#combine region and dev for easier labeling
sites.thinned<-sites.thinned %>% 
  unite(dev_region, c(region, Dev), sep="_", remove=FALSE)

#reorder the population names to match the sample order in genos
sites.thinned.sorted <- sites.thinned[order(sites.thinned$Sample_ID), ]
# Remove spaces from the sampleID column
sites.thinned.sorted$Sample_ID <- gsub(" ", "", sites.thinned.sorted$Sample_ID)

# Read sample names from the corresponding .fam file
sample_names <- read.table("genomicdata/7.urbanurchfull.new.fixed.fam", header = FALSE)[,1]  
# Assign row names to the matrix
rownames(genos) <- sample_names

bim_data <- read.table("genomicdata/7.urbanurchfull.new.bim", header = FALSE,sep = "\t")
bim_data$snp <- paste(bim_data$V1,"_",bim_data$V4, sep = "")
snp_names <- bim_data$snp


#to account for things not lining up later
if (length(snp_names) != dim(genos)[2]) {
  stop("Mismatch: The number of SNPs in the BIM file does not match the genotype matrix.")
}

######region#####
# run outflank This will run it for region
my_fst_region<- MakeDiploidFSTMat(genos, locusNames = snp_names, popNames = sites.thinned.sorted$region)

head(my_fst_region) #looks like example!
hist(my_fst_region$FST,breaks=50)
summary(my_fst_region$FST) #highest FST is higher than the mean (which is a good sign)

# estimate the SNP outliers for sites... do some trimming
outlier_region <-OutFLANK(my_fst_region,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                          NumberOfSamples =3, qthreshold = 0.01, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_region,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#how many are outliers
Pregion <- pOutlierFinderChiSqNoCorr(my_fst_region,Fstbar=outlier_region$FSTNoCorrbar,
                                dfInferred=outlier_region$dfInferred,qthreshold=0.01,Hmin=0.1)
out_region <- Pregion$OutlierFlag==TRUE #which of the SNPs are outliers?
table(out_region) #NO OUTLIERS

######tidal height#####
# run outflank
my_fst_tidal<- MakeDiploidFSTMat(genos, locusNames = snp_names, popNames = sites.thinned.sorted$Depth..m.)

head(my_fst_tidal) #looks like example!
hist(my_fst_tidal$FST,breaks=50)
summary(my_fst_tidal$FST) #highest FST is higher than the mean (which is a good sign)

# estimate the SNP outliers for sites... do some trimming
outlier_tidal <-OutFLANK(my_fst_tidal,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                          NumberOfSamples =2, qthreshold = 0.01, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_tidal,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

Ptidal <- pOutlierFinderChiSqNoCorr(my_fst_tidal,Fstbar=outlier_tidal$FSTNoCorrbar,
                                   dfInferred=outlier_tidal$dfInferred,qthreshold=0.01,Hmin=0.1)
Ptidal$q.right <- p.adjust(Ptidal$pvaluesRightTail,method="fdr")

#identify outliers, different thresholds used... pcadapt was 0.1?
outflank_outliers_tidal <- subset(Ptidal, OutlierFlag==TRUE)
outflank_tidal_names<-outflank_outliers_tidal$LocusName

out_tidal<-(Ptidal$OutlierFlag==TRUE) #which of the SNPs are outliers?
table(out_tidal) #193 outliers...?

######urb/nonurb######
#urban vs nonurban
my_fst_urb <- MakeDiploidFSTMat(genos.filt, locusNames = pcframe.filt$snp, popNames = meta.order$Dev)

outlier_urb <-OutFLANK(my_fst_urb,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                   NumberOfSamples =length(table(meta.order$Site_ID)), qthreshold = 0.01, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_urb,withOutliers=TRUE,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)
#how many are outliers
Purb <- pOutlierFinderChiSqNoCorr(my_fst_urb,Fstbar=outlier_urb$FSTNoCorrbar,
                                   dfInferred=outlier_urb$dfInferred,qthreshold=0.01,Hmin=0.1)
Purb$q.right <- p.adjust(Purb$pvaluesRightTail,method="fdr")
length(which(Purb$q.right<0.1))
out_urb <- Purb$OutlierFlag==TRUE #which of the SNPs are outliers?
table(out_urb) #165 outliers!!!

#identify outliers, different thresholds used... pcadapt was 0.1?
outflank_outliers_urb <- subset(Purb, OutlierFlag==TRUE)
write.csv(outflank_outliers_urb, "genomicdata/outflank_urb_outliers.csv")

#####look at urban/nonurban for specific regions? #####
#subset Vic, SD and LA into separate datasets then make the urban/nonurban comparison to test this... maybe later?

#############Victoria########
#subset metadata & genos matrix to only be Victoria
Vicsamples<- meta.order[meta.order$region=="Vic",]
Vicgenos <- genos[meta.order$region=="Vic",]

# run outflank
my_fst_vic <- MakeDiploidFSTMat(Vicgenos, locusNames = pcframe.filt$snp, popNames = Vicsamples$Dev)

# estimate the SNP outliers for sites... do some trimming
outlier_vic <-OutFLANK(my_fst_vic,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                          NumberOfSamples =2, qthreshold = 0.1, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_vic,withOutliers=TRUE,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#how many are outliers, qval=0.1
Pvic <- pOutlierFinderChiSqNoCorr(my_fst_vic,Fstbar=outlier_vic$FSTNoCorrbar,
                                     dfInferred=outlier_vic$dfInferred,qthreshold=0.1,Hmin=0.1)

Pvic$q.right <- p.adjust(Pvic$pvaluesRightTail,method="fdr")
length(which(Pvic$q.right<0.1)) #17 outliers

merged_data_Vic <- left_join(Pvic, bim_data, by = c("LocusName" = "snp"))
#trim short contigs
min_snps <- 400  #all others in the 1000s

merged_data_Vic <- merged_data_Vic %>%
  group_by(V1) %>%
  filter(n() >= min_snps) %>%
  ungroup()

##for all fastman plots
ovrlp<-read.csv("genomicdata/overlapping_snps_all.csv", header=TRUE)

png("Vic.manhattan_plot.png", width = 1000, height = 400, res = 150)
fastman(merged_data_Vic,chr="V1",bp="V4", p="q.right", snp="LocusName", highlight=vic_la_snps, col="blues", col2 = "greys", suggestiveline=1)
dev.off() #closes graphics device and writes file

#identify outliers, different thresholds used... 
outflank_outliers_Vic <- subset(Pvic, OutlierFlag==TRUE)
outflank_Vic_names<-outflank_outliers_Vic$LocusName

###############LA########
LAsamples<- meta.order[meta.order$region=="LA",]
LAgenos <- genos[meta.order$region=="LA",]

# run outflank
my_fst_LA <- MakeDiploidFSTMat(LAgenos, locusNames = pcframe.filt$snp, popNames = LAsamples$Dev)


# estimate the SNP outliers for sites... do some trimming
outlier_LA <-OutFLANK(my_fst_LA,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                      NumberOfSamples =2, qthreshold = 0.1, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_LA,withOutliers=TRUE,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#how many are outliers
PLA <- pOutlierFinderChiSqNoCorr(my_fst_LA,Fstbar=outlier_LA$FSTNoCorrbar,
                                 dfInferred=outlier_LA$dfInferred,qthreshold=0.1,Hmin=0.1)

PLA$q.right <- p.adjust(PLA$pvaluesRightTail,method="fdr")
length(which(PLA$q.right<0.1)) #2271!

merged_data_LA <- left_join(PLA, bim_data, by = c("LocusName" = "snp"))
#trim short contigs
min_snps <- 400  #all others in the 1000s

merged_data_LA <- merged_data_LA %>%
  group_by(V1) %>%
  filter(n() >= min_snps) %>%
  ungroup()

png("LA.manhattan_plot.png", width = 1000, height = 400, res = 150)
fastman(merged_data_LA,chr="V1",bp="V4", p="q.right", snp="LocusName", highlight=ovrlp, col="blues", col2 = "greys", suggestiveline=1)
dev.off() #closes graphics device and writes file


#identify outliers
outflank_outliers_LA <- subset(PLA, OutlierFlag==TRUE)
outflank_LA_names<-outflank_outliers_LA$LocusName

#plot(PLA$LocusName,PLA$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
#points(PLA$LocusName[out_LA],PLA$FST[out_LA],col="green")


###############SD########
#subset metadata & genos matrix 
SDsamples<- meta.order[meta.order$region=="SD",]
SDgenos <- genos[meta.order$region=="SD",]

# run outflank
my_fst_SD <- MakeDiploidFSTMat(SDgenos, locusNames = pcframe.filt$snp, popNames = SDsamples$Dev)

# estimate the SNP outliers for sites... do some trimming
outlier_SD <-OutFLANK(my_fst_SD,LeftTrimFraction=0.01,RightTrimFraction=0.01,
                       NumberOfSamples =2, qthreshold = 0.1, Hmin = 0.05)

OutFLANKResultsPlotter(outlier_SD,withOutliers=TRUE,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#how many are outliers
PSD <- pOutlierFinderChiSqNoCorr(my_fst_SD,Fstbar=outlier_SD$FSTNoCorrbar,
                                  dfInferred=outlier_SD$dfInferred,qthreshold=0.1,Hmin=0.1)

PSD$q.right <- p.adjust(PSD$pvaluesRightTail,method="fdr")
length(which(PSD$q.right<0.1)) #749

out_SD <- PSD$OutlierFlag==TRUE #which of the SNPs are outliers?
table(out_SD)

merged_data_SD <- left_join(PSD, bim_data, by = c("LocusName" = "snp"))
#trim short contigs
min_snps <- 400  #all others in the 1000s

merged_data_SD <- merged_data_SD %>%
  group_by(V1) %>%
  filter(n() >= min_snps) %>%
  ungroup()

png("SD.manhattan_plot.png", width = 1000, height = 400, res = 150)
fastman(merged_data_SD,chr="V1",bp="V4", p="q.right", snp="LocusName",highlight=la_sd_snps, col="blues", col2 = "greys",suggestiveline=1)
dev.off() #closes graphics device and writes file

#identify outliers
outflank_outliers_SD <- subset(PSD, OutlierFlag==TRUE)
outflank_SD_names<-outflank_outliers_SD$LocusName


#overlap outliers####
outframe <- data.frame(Chr=pcframe.filt$Chr,
                       Pos=pcframe.filt$Pos,
                       Snp=pcframe.filt$snp_names,
                       my_fst_urb=Purb$FST[match(pcframe.filt$snp_names,Purb$LocusName)],
                       OutflankUrb=Purb$q.right[match(pcframe.filt$snp_names,Purb$LocusName)],
                       my_fst_SD=PSD$FST[match(pcframe.filt$snp_names,PSD$LocusName)],
                       SD.q=PSD$q.right[match(pcframe.filt$snp_names,PSD$LocusName)],
                       my_fst_LA=PLA$FST[match(pcframe.filt$snp_names,PLA$LocusName)],
                       LA.q=PLA$q.right[match(pcframe.filt$snp_names,PLA$LocusName)],
                       my_fst_vic=Pvic$FST[match(pcframe.filt$snp_names,Pvic$LocusName)],
                       Vic.q=Pvic$q.right[match(pcframe.filt$snp_names,Pvic$LocusName)])

OutflankUrb <- na.omit(data.frame(pcframe.filt[,4:10],outframe[,c(5,7,9,11)]))
dim(OutflankUrb)

All.tf <- OutflankUrb<0.1

#all pcadapt vs outflank outliers
plot(euler(All.tf[,c("PC.full","OutflankUrb")],shape="ellipse"),quantities=list(cex = 1.5),
     edges="black")

#region specific outliers with outflank + pcadapt-- select columns to include in
plot(euler(All.tf[,c(1,9:11)],shape="ellipse"),quantities=list(cex = 1.5),
     edges="black")

#only regions
#region specific outliers with outflank + pcadapt-- select columns to include in
VeD<-
  plot(euler(All.tf[,c(9:11)],shape="ellipse"),quantities=list(cex = 4),
     edges="black")
ggsave("venndiagram.png", VeD, width=30, height=20, units = "cm")



#top 5% SNPs for regions and methods
out.quant <- data.frame(All=outframe$my_fst_urb>quantile(outframe$my_fst_urb,0.95,na.rm=T),
                        SD=outframe$my_fst_SD>quantile(outframe$my_fst_SD,0.95,na.rm=T),
                        LA=outframe$my_fst_LA>quantile(outframe$my_fst_LA,0.95,na.rm=T),
                        Vic=outframe$my_fst_vic>quantile(outframe$my_fst_vic,0.95,na.rm=T))
plot(euler(na.omit(out.quant),shape="ellipse"),quantities=list(cex = 1.5),
     edges="black")

##upset plot####
library(UpSetR)
#upset plot for top 5% SNPs
expressionInput <- c("All"= 4101,"Vic" = 8292, "LA" = 6409, "SD" = 8525, "Vic&LA" = 0, "Vic&SD" = 382,
                     "LA&SD" = 323, "All&Vic"=1506, "All&SD"=1241, "All&SD&Vic"=170, "All&LA&SD"=240, "LA&All"=3386)

upset(fromExpression(expressionInput), order.by="degree",
      matrix.color="red")

#upset plot for outlier overlap
expressionInput1 <- c("pcadapt"=44, "Vic" = 16, "LA" = 1569, "SD" = 417, "Vic&LA" = 1, "Vic&SD" = 0,
                     "LA&SD" = 5, "pcadapt&LA"=4)

upset(fromExpression(expressionInput1), order.by="degree",
      matrix.color="purple")


###no longer will work I dont think??
# Intersection of Victoria & LA
vic_la_snps <- intersect(outflank_Vic_names, outflank_LA_names)
# Intersection of Victoria & San Diego--none
# Intersection of LA & San Diego
la_sd_snps <- intersect(outflank_LA_names, outflank_SD_names)
#pcadapt x LA intersection is that even needed?

# Create a data frame with columns of equal length
overlapping_snps <- data.frame(vic_la_snps, la_sd_snps)

# Save to CSV
write.csv(overlapping_snps, "genomicdata/overlapping_snps.csv", row.names = FALSE)



###Polygenic scores- RB####
phen <- as.numeric(as.factor(meta.order$Dev))
##Read imputed genotypes
imp.g <- read.table("genomicdata/7.urbanurchfull.new.lfmm_imputed.lfmm",sep=" ")

###Run LFMM on all samples
lfmm <- lfmm2(input=as.matrix(imp.g),
              env=as.matrix(phen),
              K=2,
              effect.sizes=T)
pv <- lfmm2.test(object=lfmm,
                 input=as.matrix(imp.g),
                 env=as.matrix(phen),
                 linear=T,genomic.control=T)
plot(-log10(pv$pvalues))
q.vals <- p.adjust(pv$pvalues,method="fdr")
length(which(q.vals<0.1)) #no outliers

lfmm.frame <- data.frame(sites, pvalue=pv$pvalues,qvalue=q.vals)
write.table(lfmm.frame,"ubanurchin.lfmm.results.txt",quote=F,row.names=F)
n=100
res.frame <- data.frame(LA.n=NA,
                        LA.u=NA,
                        SD.n=NA,
                        SD.u=NA,
                        VIC.n=NA,
                        VIC.u=NA,
                        p=NA)

##subset training set currently 60% (110 individuals)
for (i in 1:n) {
  print(i)
  train <- as.matrix(sample(1:nrow(meta.order),110))
  test <- setdiff(1:nrow(meta.order),train)
  train.g <- imp.g[train,]
  train.e <- phen[train]
  valid.g <- imp.g[test,]
  valid.e <- phen[test]
  
  ##Run LFMM
  lfmm <- lfmm2(input=as.matrix(train.g),
                env=as.matrix(train.e),
                K=2,
                effect.sizes=T)
  
  
  ## Effect sizes
  b.values <- lfmm@B[,1]
  
  #Prediction
  pred <- scale(valid.g,scale=F) %*% matrix(b.values)
  
  frame <- data.frame(Dev=meta.order$Dev[test],
                      Reg=meta.order$region[test],
                      Score=pred)
  
  agg <- aggregate(frame$Score,list(frame$Dev,frame$Reg),mean,na.rm=T)
  res.frame[i,1:6] <- agg$x
  mod <- anova(lm(Score~Dev+Reg, data=frame))
  res.frame[i,7] <- mod$`Pr(>F)`[1]
}

gath <- gather(res.frame,"metric","value")
write.csv(gath,"Polygenic_results.csv")

#How many bootstraps are significant?
length(which(res.frame$p<0.05))

#Make boxplot
gath.format <- gath %>% filter(metric!="p") %>% 
  separate(metric,into=c("Region","Dev"),sep="\\.") %>%
  mutate(Dev=case_when(Dev=="n"~"nonurban",
                       Dev=="u"~"urban"))
ggplot(gath.format,aes(x=Region,y=value,fill=Dev)) + 
  geom_boxplot() + theme_bw() +
  ylab("Polygenic Score") + scale_fill_viridis(discrete=T,begin=0.3,end=0.7)



###Null for polygenic scores using randomized environment 

n=100
res.frame <- data.frame(LA.n=NA,
                        LA.u=NA,
                        SD.n=NA,
                        SD.u=NA,
                        VIC.n=NA,
                        VIC.u=NA,
                        p=NA)

##subset training set currently 60% (110 individuals)
for (i in 1:n) {
  print(i)
  train <- as.matrix(sample(1:nrow(meta.order),110))
  test <- setdiff(1:nrow(meta.order),train)
  train2 <- sample(1:nrow(meta.order),110) #for randomizations
  test2 <- setdiff(1:nrow(meta.order),train2) #for randomizations
  train.g <- imp.g[train,]
  train.e <- phen[train2] #for randomizations
  valid.g <- imp.g[test,]
  valid.e <- phen[test2]
  
  ##Run LFMM
  lfmm <- lfmm2(input=as.matrix(train.g),
                env=as.matrix(train.e),
                K=2,
                effect.sizes=T)
  
  
  ## Effect sizes
  b.values <- lfmm@B[,1]
  
  ## Prediction
  pred <- scale(valid.g,scale=F) %*% matrix(b.values)
  
  frame <- data.frame(Dev=meta.order$Dev[test],
                      Reg=meta.order$region[test],
                      Score=pred)
  
  agg <- aggregate(frame$Score,list(frame$Dev,frame$Reg),mean,na.rm=T)
  res.frame[i,1:6] <- agg$x
  mod <- anova(lm(Score~Dev+Reg, data=frame))
  res.frame[i,7] <- mod$`Pr(>F)`[1]
}

#How many bootstraps are significant?
length(which(res.frame$p<0.05))
gath <- gather(res.frame,"metric","value")
write.csv(gath,"Polygenic_randomizations.csv")



#Make polygenic boxplot ####

##if data was totally random and not driven by urban/nonurban
gath<- read.csv("genomicdata/Polygenic_randomizations.csv")
gath.format <- gath %>% filter(metric!="p") %>% 
  separate(metric,into=c("Region","Dev"),sep="\\.") %>%
  mutate(Dev=case_when(Dev=="n"~"nonurban",
                       Dev=="u"~"urban"))
#p1<-
  ggplot(gath.format,aes(x=Region,y=value,fill=Dev)) + 
  geom_boxplot() + theme_box() +
  ylab("Polygenic Score") +ggtitle("Random Grouping") + scale_fill_manual(values = c("nonurban"="orchid2" ,"urban" = "tan"))
#ggsave("polygenic_random.png", p1, width=25, height=20, units = "cm")

#driven by urban/nonurban
nonrandom<- read.csv("genomicdata/Polygenic_results.csv")
nonrandom.format <- nonrandom %>% filter(metric!="p") %>% 
  separate(metric,into=c("Region","Dev"),sep="\\.") %>%
  mutate(Dev=case_when(Dev=="n"~"nonurban",
                       Dev=="u"~"urban"))
p2<-
ggplot(nonrandom.format,aes(x=Region,y=value,fill=Dev)) + 
  geom_boxplot() + theme_box() +
  ylab("Polygenic Score") +ggtitle("D")+ scale_fill_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))
ggsave("polygenic_nonrandom.png", p2, width=20, height=15, units = "cm")
