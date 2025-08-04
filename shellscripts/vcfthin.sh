#!/bin/bash

#SBATCH --job-name=vcfthin
#SBATCH --mem=100G
#SBATCH -o /group/rbaygrp/armstrong/urbanurchins/logs/out-%A.%a_vcfthin.txt
#SBATCH -e /group/rbaygrp/armstrong/urbanurchins/logs/error-%A.%a_vcfthin.txt
#SBATCH --time=2-2:00:00
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -p high2

#to run sbatch vcfthin.sh
###
#load modules and set working directory
module load vcftools

#First, remove low coverage individuals (vcftools), code below

##Then, remove bad SNPs
###Lenient-est filter
vcftools --gzvcf /group/rbaygrp/armstrong/urbanurchins/vcf_results/ubanurchinfull_raw.vcf.gz \
--remove /group/rbaygrp/armstrong/urbanurchins/vcf_results/low.cov.samps.txt \
--min-alleles 2 --max-alleles 2 \
--remove-indels \
--recode \
--minGQ 20 \
--max-missing 0.60 --maf 0.01 --out /group/rbaygrp/armstrong/urbanurchins/vcf_results/3.ubanurchinfull_prelimfilter60

##Thin data
vcftools --vcf /group/rbaygrp/armstrong/urbanurchins/vcf_results/3.ubanurchinfull_prelimfilter60.recode.vcf \
--thin 25000 \
--site-mean-depth \
--out /group/rbaygrp/armstrong/urbanurchins/vcf_results/4.ubanurchinfull_prelimfilter60_thinned
#Assess depth to adjust parameters below
##Plot histogram - see the tutorial I sent


