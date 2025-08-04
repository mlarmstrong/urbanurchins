#!/bin/bash

#SBATCH --job-name=vcfthin2
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH -o /group/rbaygrp/armstrong/urbanurchins/logs/out-%A.%a_vcfthin2.txt
#SBATCH -e /group/rbaygrp/armstrong/urbanurchins/logs/error-%A.%a_vcfthin2.txt
#SBATCH --time=1-12:00:00
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -p high2

#Take what you learned from the histogram and apply here!

#summary(var_depth$mean_depth)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4.770   5.984   6.492   6.912   7.279 106.230
#min is usually good to put at 10 but going to do 8, max depth usually avg*2 so ~20 to be flexible

vcftools --vcf /group/rbaygrp/armstrong/urbanurchins/vcf_results/3.ubanurchinfull_prelimfilter10.recode.vcf \
--recode \
--maf 0.01 --min-meanDP 8 --max-meanDP 20 \
--max-missing 0.70 \
--out /group/rbaygrp/armstrong/urbanurchins/vcf_results/5.ubanurchinfull_fullyfiltered10

vcftools --vcf /group/rbaygrp/armstrong/urbanurchins/vcf_results/5.ubanurchinfull_fullyfiltered10.recode.vcf \
--thin 25000 \
--recode \
--out /group/rbaygrp/armstrong/urbanurchins/vcf_results/5.ubanurchinfull_fullyfiltered10_thin
