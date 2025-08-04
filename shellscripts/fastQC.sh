#! /bin/bash
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=QC
#SBATCH --mem=60G
#SBATCH --partition=med
#SBATCH --time=3-1:00:00
set -e
set -x
# To run script type:  sbatch fastQC.sh

## Change directories to where the fastq files are located
cd /group/rbaygrp/armstrong/urbanurchins/usemedata
# Call fastp package
module load fastqc

#RUN FOR IT MARTY
fastqc -t 5 *.fq.gz -o /group/rbaygrp/armstrong/urbanurchins/qual
