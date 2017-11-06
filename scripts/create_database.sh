#!/bin/bash
set -e

## Download Phage SEED
mkdir phage_seed_jun2017
cd phage_seed_jun2017
wget "http://www.phantome.org/Downloads/proteins/all_sequences/phage_proteins_1496314803.fasta.gz"
gunzip phage_proteins_1496314803.fasta.gz
mv phage_proteins_1496314803.fasta PHAGE_SEED_2017-06-01

## Make the blast database
makeblastdb -in PHAGE_SEED_2017-06-01 -dbtype prot

## Download the taxonomy information
mkdir tax
cd tax
wget "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
tar -xvf taxdump.tar.gz
rm taxdump.tar.gz
cd ../
