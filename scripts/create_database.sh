#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p BabyViromeDB
cd BabyViromeDB

## Download SEED subsystem information
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "ftp://ftp.theseed.org/subsystems/subsystems2peg.gz"
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "ftp://ftp.theseed.org/subsystems/subsys.txt"
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "ftp://ftp.theseed.org/genomes/SEED/all.faa.gz"

## Download Phage SEED
wget "http://www.phantome.org/Downloads/proteins/all_sequences/phage_proteins_1496314803.fasta.gz"

## Set everything up
gunzip phage_proteins_1496314803.fasta.gz
mv phage_proteins_1496314803.fasta Phage_SEED_2017-06-01
gunzip subsystems2peg.gz
gunzip all.faa.gz
mv all.faa SEED_original
"${SCRIPT_DIR}/rebuild_seed_headers.pl" -s2 subsystems2peg -s subsys.txt -f SEED_original
mv SEED_original_updated.fasta SEED
rm SEED_original

## Make the blast databases
makeblastdb -in SEED -dbtype prot
makeblastdb -in Phage_SEED_2017-06-01 -dbtype prot

## Download the taxonomy information
mkdir tax
cd tax
wget "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
tar -xvf taxdump.tar.gz
rm taxdump.tar.gz
mv nodes.dmp ../
mv names.dmp ../
cd ../
rm -rf tax
cd ../
