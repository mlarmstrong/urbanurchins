##pop structure & pca urchin script
#Madison Armstrong
#12/15/24
#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "white", fill = NA, size = 1)) 
}
##SET UP FOR DIFF ANALYSES ####

#tess3r packages
#devtools::install_github("bcm-uga/TESS3_encho_sen")
library(tess3r)
library(vcfR)
library(maps)
library(tidyverse)
# install pophelper package from GitHub
#remotes::install_github('royfrancis/pophelper')
library(pophelper)
library(gridExtra)
library(ggplot2)

#pca
library(devtools)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(vcfR)
library(viridis)
library(dplyr)
library(dichromat)
library(tidyverse)

#set wd
setwd("~/Desktop/urbanurchins")
########


#####tess3r####
#read vcf file
vcf <- read.vcfR("genomicdata/5.ubanurchinfull_fullyfiltered10.new_thin.recode.vcf")
genos_raw<-extract.gt(vcf) #extract genotypes, format is 0/0, 0|1, 0/1 etc

# Replace genotype codes with numeric values
genos_col <- as.matrix(apply(genos_raw, c(1, 2), function(x) {
  if (x %in% c("0/0", "0|0")) return(0)    # Homozygous reference → 0
  if (x %in% c("0/1", "1/0", "0|1", "1|0")) return(1)  # Heterozygous → 1
  if (x %in% c("1/1", "1|1")) return(2)    # Homozygous alternate → 2
  return(NA)             # Handle missing or other cases
}))
# Check the cleaned matrix
#View(genos_col) #Genotypes are encoded as 0,1,2 for diploids
#need to flip rows and columns
dim(genos_col) #indiv should be on left, genotypes on right

# Flip the genotype matrix (transpose)
genos<- t(genos_col) 
dim(genos) #now it is right!! 183, 19299
any(is.na(genos)) #TRUE so there are NA values in the dataset

#and file with coords
sites.thinned<-read.csv('metadata/urbanurchins_metadata_thinned.csv', header=TRUE, sep=',')
coords<-as.matrix(cbind(sites.thinned$Latitude, sites.thinned$Longitude))
sample_id<-as.matrix(cbind(sites.thinned$Sample_ID, sites.thinned$Latitude, sites.thinned$Longitude))

#estimating ancestry coeffs
#The X argument refers to the genotype matrix
#the coord argument corresponds to the geographic coordinates
#K is the number of clusters or ancestral population
#and openMP.core.num is the number of processes used by the multi-threaded program. or do rep=some number
#use "try" to print errors
tess3.obj <- tess3(X=genos, coord =as.matrix(sites.thinned[,c("Longitude", "Latitude")]), K=1:8,
                   ploidy=2, lambda=0, keep="best",rep=10)
#print(tess3.obj)
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score") 
#best K= when curve starts to flatten, but can also try to estimate statistically

# retrieve tess3 Q matrix for K = 2 clusters 
q.matrix <- qmatrix(tess3.obj, K = 2)

# STRUCTURE-like barplot for the Q-matrix 
#combine region and dev for easier labeling
sites.thinned<-sites.thinned %>% 
  unite(dev_region, c(region, Dev), sep="_", remove=TRUE)

barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = sites.thinned$dev_region,las = 3, cex.axis = .4) 
#labeled with sample_IDs but need more organization!

#are individuals grouping separately based on the qmatrix? Test as "pop" for the color

# Install and load adegenet
#install.packages("adegenet")
library(adegenet)

# Convert Q-matrix to a genind object
genind_obj <- df2genind(q.matrix, pop = as.factor(sites.thinned$Site_ID))

###look at individual regions####
sites.thinned<-read.csv('metadata/urbanurchins_metadata_thinned.csv', header=TRUE, sep=',')
#Victoria
sites.Vic<- sites.thinned[sites.thinned$region=="Vic",]
genos.Vic<- genos[sites.thinned$region=="Vic",]
tess3.vic<- tess3(X=genos.Vic, coord =as.matrix(sites.Vic[,c("Longitude", "Latitude")]), K=1:8,
                   ploidy=2, lambda=0, keep="all",rep=10)
#print(tess3.obj)
#plot(tess3.vic, pch = 19, col = "purple",
   #  xlab = "Number of ancestral populations",
    # ylab = "Cross-validation score")
q.matrix.v <- qmatrix(tess3.vic, K = 2, rep=2)


barplot(q.matrix.v, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix.v), labels = sites.Vic$Site_ID,las = 3, cex.axis = 1) 

#LA
sites.LA<- sites.thinned[sites.thinned$region=="LA",]
genos.LA<- genos[sites.thinned$region=="LA",]
tess3.LA<- tess3(X=genos.LA, coord =as.matrix(sites.LA[,c("Longitude", "Latitude")]), K=1:5,
                  ploidy=2, lambda=0, keep="best",rep=10)
#print(tess3.obj)
plot(tess3.LA, pch = 19, col = "pink",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score") 
q.matrix.la <- qmatrix(tess3.LA, K = 4)

barplot(q.matrix.la, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix.la), labels = sites.LA$Site_ID,las = 3, cex.axis = 1) 


#SD
sites.SD<- sites.thinned[sites.thinned$region=="SD",]
genos.SD<- genos[sites.thinned$region=="SD",]
tess3.SD<- tess3(X=genos.SD, coord =as.matrix(sites.SD[,c("Longitude", "Latitude")]), K=1:5,
                 ploidy=2, lambda=0, keep="best",rep=10)

plot(tess3.SD, pch = 19, col = "orange",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score") 
q.matrix.SD <- qmatrix(tess3.SD, K = 4)

barplot(q.matrix.SD, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix.SD), labels = sites.SD$Site_ID,las = 3, cex.axis = 1) 

##pophelper####
# Convert Q-matrix to a genind object
genind_obj <- df2genind(q.matrix, pop = as.factor(sites.thinned$Site_ID))

#use pophelper! https://www.royfrancis.com/pophelper/articles/index.html#plotting-1

# convert TESS 3 object to qlist for pophelper
qlist <- readQTess3(t3list=tess3.obj)
is.qlist(qlist=qlist)#verify format

#unique(sites.thinned$Site_ID) #get list of sites
poporder=c("Vic_urban","Vic_nonurban","LA_urban","LA_nonurban","SD_urban", "SD_nonurban")
labs <- c("K=2", "K=3")
depth=c("intertidal", "subtidal")


#combine region and dev for easier labeling
labset<-sites.thinned

labset<-labset[,3:8] #pull data that we want from sites.thinned

slist1 <- alignK(qlist[2:3],type="auto")
verifyGrplab(grplab=labset[,c("dev_region", "Site_ID")])#verify format

#plot samples separated by region_dev
p1<-plotQ(slist1,imgoutput="join",returnplot=TRUE,exportplot=FALSE,basesize=11,
            splabsize=7,height=7,
            grplab=labset[,c("dev_region", "Site_ID")],subsetgrp=poporder,
            grplabsize=3,linesize=1,pointsize=3,splab=labs,grplabangle=0,
            grplabheight = 5)

#add intertidal/subtidal labels
p2<-plotQ(slist1,imgoutput="join",returnplot=T,exportplot=F,basesize=11,
          splabsize=7,height=7,
          grplab=labset[,c("dev_region", "Depth..m.")],subsetgrp=poporder,
          grplabsize=4,linesize=1,pointsize=3,splab=labs,grplabangle=0,
          grplabheight = 5) 
grid.arrange(p1$plot[[1]], p2$plot[[1]])
popstruc<-grid.arrange(p1$plot[[1]])

ggsave("popstructure.png", popstruc, width=30, height=20, units = "cm")

###PCA #####
# load data and convert to gds format

##SEQVCF2GDS large file
##convert vcf to gds and open file: using seqVCF2GDS
#better for large data and combining vcf files
vcf.fn <- "genomicdata/5.ubanurchinfull_fullyfiltered10.new_thin.recode.vcf"
seqVCF2GDS(vcf.fn, "ubanurchinfull.gds")
genofileseq <- seqOpen("ubanurchinfull.gds")

pca<-snpgdsPCA(genofileseq, autosome.only=FALSE)


names(pca) #look at column headers
#"sample.id" "snp.id"    "eigenval"  "eigenvect" "varprop"   "TraceXTX"  "Bayesian"  "genmat"   

  
#variance explained by each PC
pc.percent<-pca$varprop*100
plot(pc.percent)
print(pc.percent) #to get the values for the pca plot axes

#look at pc.percent values in environment

##Plotting PCAs
tab<- data.frame(sample.id=pca$sample.id,
                 EV1=pca$eigenvect[,1],
                 EV2=pca$eigenvect[,2],
                 EV3=pca$eigenvect[,3],
                 EV4=pca$eigenvect[,4])


#Population meta
sites.thinned<-read.csv('metadata/urbanurchins_metadata_thinned.csv', header=TRUE, sep=',')

#arrange by latitude
sites.thinned<-sites.thinned %>% 
  arrange(Latitude) %>%
  unite(dev_region, c(region, Dev), sep="_", remove=FALSE)
sites.thinned$Site_ID <- factor(sites.thinned$Site_ID, levels = unique(sites.thinned$Site_ID))


#PLOTS
#ensure shapes separate urban/nonurban
shapes<-ifelse(sites.thinned$dev_region)
ggplot(tab, aes(x = EV2, y = EV1,fill=as.factor(sites.thinned$Dev))) +geom_point(size = 5, stroke = 1, aes(color = I("black")))
#pca<-
ggplot(tab, aes(x = EV2, y = EV1, fill=as.factor(sites.thinned$Dev), shape=sites.thinned$region))+
  geom_point(size = 5, stroke = 1, aes(color = I("black"))) +
  labs(x = 'PC1 0.7256147%', y='PC2 0.7241320%', fill = "Development",  # Legend title for fill
       shape = "Region")+
  theme_minimal(base_size = 18)+
  scale_fill_manual(name="Development",values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+
  scale_shape_manual(name="Region",values = c(21, 22, 24)) +
  guides(color = "none")+
  theme_box()

ggsave("PCA_neutral.png", pca, width=30, height=20, units = "cm")

ggplot(tab, aes(x = EV2, y = EV1, shape=sites.thinned$dev_region))+
  geom_point(size = 5)+
  labs(x = 'PC1 0.7256147%', y='PC2 0.7241320%')+
  theme_minimal(base_size = 18)+
  theme(panel.background = element_rect(fill = "white"))

ggplot(tab, aes(x = EV4, y = EV3, color=as.factor(sites.thinned$Dev), shape=sites.thinned$dev_region))+
  geom_point(size = 5)+
  labs(x = 'PC3 0.7175024%', y='PC4 0.7039728%')+
  theme_minimal(base_size = 18)+
  theme(panel.background = element_rect(fill = "white"))

