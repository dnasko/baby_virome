#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --job-name=Example
#SBATCH --qos=large
#SBATCH --partition=dpart
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20000

## HINTS!!
## time is HH:MM:SS
## give your job an actual name, bro
## qos/partition should be fine for umiacs, change if not UMIACS
## ntasks-per-node is the number of CPUs you need. It's probably one unless
##                 you have a multi-threaded program
## mem is memory in MB (I know, dumb). 20000 is 20 GB

hostname; date;
echo " Machine has $SLURM_CPUS_ON_NODE CPU cores"

## Put whatever command you want down here...


date
