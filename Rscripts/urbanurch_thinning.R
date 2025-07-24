#popgen thinning urchin script
#Madison Armstrong
#12/13/24
#first histogram in order to know how to best thin the vcf data

library(tidyverse)

setwd("~/Desktop/urbanurchins")
#https://speciationgenomics.github.io/filtering_vcfs/ 
#start at tutorial "variant mean depth"

#GQ=10, max-missing=0.7
var_depth <- read_delim("genomicdata/4.ubanurchinfull_prelimfilter10_thinned.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)


a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() #looks similar to tutorial

summary(var_depth$mean_depth)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.770   5.984   6.492   6.912   7.279 106.230  

#min is usually good to put at 10, max depth usually avg*2 so ~14 ?
#redraw to ignore outliers 
a + theme_light() + xlim(0, 100)

###OLD#####
#GQ=20, max-missing=0.7
var_depth <- read_delim("genomicdata/4.ubanurchinfull_prelimfilter_thinned.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)


a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() #looks similar to tutorial

summary(var_depth$mean_depth)
#   Min. 1st Qu.  Median  Mean    3rd Qu.    Max. 
# 7.525  10.847  11.661  12.517  12.940   197.372 

#min is usually good to put at 10, max depth usually avg*2 so ~26
#redraw to ignore outliers 
a + theme_light() + xlim(0, 100)

