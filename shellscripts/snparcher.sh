#!/bin/bash

#SBATCH --job-name=sm
#SBATCH --mem=100G
#SBATCH --ntasks=8
#SBATCH -o /group/rbaygrp/armstrong/urbanurchins/logs/out-%A.%a_snakemake.txt
#SBATCH -e /group/rbaygrp/armstrong/urbanurchins/logs/error-%A.%a_snakemake.txt
#SBATCH --time=10-12:00:00
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -p high2

module load conda/base/latest
source activate snakemake --unlock
mamba activate snparcher

cd /group/rbaygrp/armstrong/urbanurchins/snpArcher
snakemake --snakefile workflow/Snakefile --workflow-profile ./profiles/slurm

