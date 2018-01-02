#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00
#SBATCH --job-name=MGAssembly
#SBATCH --qos=large
#SBATCH --partition=dpart
#SBATCH --ntasks-per-node=24
#SBATCH --mem=120000

hostname; date;
echo " Machine has $SLURM_CPUS_ON_NODE CPU cores"

ROOT="/cbcbhomes/jchopyk/projects/conserve_year_1"

## This is the script down here.
## So long as your reads are in the /reads/ folder and each file ends in "_R#.fastq", then
##  you only need to update the SAMPLE variable below:
SAMPLE="RWW_july"

## Command
time $ROOT/metagenome_assembly_pipeline.sh \
    $ROOT/reads/${SAMPLE}_R1.fastq \
    $ROOT/reads/${SAMPLE}_R2.fastq \
    $SAMPLE \
    $ROOT/assemblies \
    24
date
