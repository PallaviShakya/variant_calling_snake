#! /bin/bash -l
#SBATCH -J F2s_Variant_calling
#SBATCH -e out.erri
#SBATCH -o Test-Job-%j.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pshakya@ucdavis.edu
#SBATCH --account=datalabgrp
#SBATCH --partition=med2
#SBATCH --time=20-7:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB

#Make things fail on errors

set -o errexit
set -x

conda activate variant_calling

# snakemake -n # Dry run first

snakemake --unlock --rerun-incomplete && snakemake --cores all 

#extra stuff used: 
#samtools faidx to index reference
#

#gatk CreateSequenceDictionary -R ref.fasta ## create dict file


