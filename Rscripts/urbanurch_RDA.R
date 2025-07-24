#urbanurch_RDA
#Madison Armstrong
#3/31/25
#tutorials: https://popgen.nescent.org/2018-03-27_RDA_GEA.html, https://marineomics.github.io/RDAtraitPredictionTutorial.html#Load_the_data

setwd("~Desktop/urbanurchins")
library(psych) # Used to investigate correlations among predictors
library(vegan) # Used to run RDA
library(dplyr) # Used for data cleaning
library(tidyr) # Used for data cleaning
library(ggplot2) # Used for plotting
library(plotly) # Used with ggplot
library(ggpubr) # Used with ggplot
library(sf) # Used to make maps
library(rnaturalearth) # Used to make maps
library(rnaturalearthdata) # Used to make maps
library(maps) # Used to make maps
library(ggrepel) # Used to make maps
library(ggspatial) # Used to make maps
library(LEA) #need to install with Biomanager
library(vegan)

###read in data####
#RDA requires no missing data/complete data frames, so need to impute using most common genotype at each SNP across all individuals
#trim by MAF=0.05

#read in imputed matrix from LFMM, 7.urbanurchfull.new.lfmm_imputed.lfmm
gen.imp<- read.lfmm("genomicdata/7.urbanurchfull.new.lfmm_imputed.lfmm")

#read in metadata
sites.thinned<-read.csv('metadata/urbanurchins_metadata_thinned.csv', header=TRUE, sep=',')
#combine region and dev for easier labeling
sites.thinned<-sites.thinned %>% 
  unite(dev_region, c(region, Dev), sep="_", remove=FALSE)

#reorder the population names to match the sample order in genos
sites.thinned.sorted <- sites.thinned[order(sites.thinned$Sample_ID), ]
# Remove spaces from the sampleID column
sites.thinned.sorted$Sample_ID <- gsub(" ", "", sites.thinned.sorted$Sample_ID)
sites.test<-
  sites.thinned.sorted %>% select(Sample_ID, Site_ID, Latitude, Longitude, region, Dev, Depth..m.)

sites.test$Site_ID <- as.factor(sites.test$Site_ID)
sites.test$Dev <-as.factor(sites.test$Dev)
sites.test$Depth..m. <-as.factor(sites.test$Depth..m.)

#look at correlations
pairs.panels(sites.test, scale=T)

sites.final<-
  sites.test %>% select(Latitude, Dev, Depth..m.)

#look at correlations
pairs.panels(sites.final, scale=T)

#RUN RDA ####
urch.rda <- rda(gen.imp ~ Latitude + Dev + Depth..m., data=sites.test)
urch.rda

RsquareAdj(urch.rda)
summary(eigenvals(urch.rda, model = "constrained"))

# MODEL SELECTION####
# first build empty model for use in ordistep
urch.rda.0 <- rda(gen.imp ~ 1, sites.test)
urch.rda.0

# backwards model selection with ordistep
ordistep(urch.rda, urch.rda.0, direction="back", permutations = how(nperm=1000), step = 100)

###run models#
##need to run on pompeii worm, I don't have enough memory 
signif.full <- anova.cca(urch.rda, parallel=getOption("mc.cores"), permutations = how(nperm=200))
signif.full

signif.axis <- anova.cca(urch.rda, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=200))
signif.axis #RDA 1 and 2 

vif.cca(urch.rda)

#run model without Lat
urch.rda.L <- rda(gen.imp ~ Dev + Depth..m., data=sites.test)
urch.rda.L

#compare the two models --pompeii worm
#cmpr<-anova(urch.rda.L, urch.rda, permutations = 999)
###including latitude doesn’t significantly explain more of the variation##

#signif.L <- anova.cca(urch.rda.L, parallel=getOption("mc.cores"), permutations = how(nperm=200))
#signif.L #NO LAT IS ALSO SIG

#signif.axis.L <- anova.cca(urch.rda.L, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=200))
#signif.axis.L

#sig.term.L <-anova.cca(urch.rda.L, by="term", parallel=getOption("mc.cores"), permutations = how(nperm=200))
#sig.term.L

# Visualizing site scores
# loadings by site
loadings.site <- scores(urch.rda.L, choices =c(1:3), display = 'sites')

# get these into the same df
# first set row names as first column
loadings.sites <- cbind(rownames(loadings.site), data.frame(loadings.site, row.names=NULL))
colnames(loadings.sites) <- c('name','RDA1','RDA2','RDA3')
dim(loadings.sites)

loadings.sites.combined <- cbind(loadings.sites, sites.thinned.sorted)
dim(loadings.sites.combined) 
write.csv(loadings.sites.combined, "loadings.sites.csv")

##base plotting
plot(urch.rda.L, scaling=3) 
#scaling 2 shows effects of explanatory variables
#scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues
text(urch.rda.L, display="bp", scaling=3) 

###CCGP RDA plotting ####
#color by dev, shape by region -- green vs grey; light blue vs brown (urban=gross)
#this is how my RDA should look
ggplot(sites.complete,aes(x=RDA1,y=RDA2,color=Dev)) + geom_point(aes(shape=region)) + theme_bw() 

# Extract scores (scaling = 2 for interpretation)
site_scores <- scores(urch.rda.L, display = "sites", scaling = 2) #pull scores
rownames(site_scores) <- sites.thinned.sorted$Sample_ID #name with sample_ID
site_scores <- as.data.frame(site_scores) %>% 
  mutate(Site = rownames(.))

#extract SNP scores
snps_scores <- scores(urch.rda.L, display = "species",choices = 1:2)
snps_scores <- as.data.frame(snps_scores)
snps_scores$SNP <- rownames(snps_scores)


#arrows
env_arrows_raw <- scores(urch.rda.L, display = "bp", scaling = 2)
# Rescale for plotting
mult <- vegan:::ordiArrowMul(env_arrows_raw)
env_arrows <- as.data.frame(env_arrows_raw * mult)

# Add variable names for plotting
env_arrows$Variable <- rownames(env_arrows)

#metadata
site_scores <- site_scores %>%
  left_join(sites.thinned.sorted, by = c("Site" = "Sample_ID"))

##plot RDA!

#"nonurban"="#99d1ec" ,"urban" = "#665d4b"

# arrows showing up short for some reason?
# 1. Compute correlation strength (arrow length)
env_arrows$cor_strength <- sqrt(env_arrows$RDA1^2 + env_arrows$RDA2^2)
base_scale <- 2 

# 3. Multiply each vector by its own correlation strength * base_scale
env_arrows$RDA1_scaled <- env_arrows$RDA1 * env_arrows$cor_strength * base_scale
env_arrows$RDA2_scaled <- env_arrows$RDA2 * env_arrows$cor_strength * base_scale

##more informative labels
env_arrows$Label <- c("Urban", "Depth")  # Example custom labels

ubanRDA<-
  ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, shape=region, fill = Dev),
             size = 3) +
    geom_jitter(data = snps_scores, aes(x = RDA1, y = RDA2), color = "black", width = 0.05, height = 0.05, alpha = 0.3, size = 2, shape=1)+
    geom_segment(data = env_arrows, 
               aes(x = 0, y = 0, xend = RDA1_scaled, yend = RDA2_scaled),
               arrow = arrow(length = unit(0.5, "cm")), 
               color = "black", linewidth = 0.8) +
  geom_text(data = env_arrows, 
            aes(x = RDA1_scaled* 1.2, y = RDA2_scaled* 1.2, label = Label), 
            color = "black", size = 5)+
  scale_shape_manual(values=c("LA"=21, "SD"=22, "Vic"=24))+
  scale_fill_manual(values=c("urban"="#665d4b", "nonurban"="#99d1ec")) +
  theme_box() +
  coord_equal() +
  labs(x = "RDA1", y = "RDA2", title = "RDA") +
  theme(legend.position = "right")


ggsave("RDA.urb.png", ubanRDA, width=20, height=15, units = "cm")

